# aq.py

# Solve 2-D aquifer pumping optimization

# Author: Michael Skarbek (maskarb@gmail.com)
# Version: 0.1.0

import itertools
import math
import random
import gurobipy as gu
import matplotlib.pyplot as plt
import numpy as np

from .report import lindo

tupledict = gu.tupledict  # pylint: disable=E1101


def _make_zero_tuple_dict(x, y):
    return tupledict((tup, 0) for tup in itertools.product(range(1, x), range(1, y)))


class Grid(object):
    """
    The Grid class represents the 2-dimensional aquifer. The class is
    instantiated with width and height and the remaining attributes are
    added with adder methods. The order in which attributes are added
    is important. Check the method documentation if errors are raised.

    Attributes:
        width (int): Number of grid cells width
        height (int): Number of grid cells height

    """
    def __init__(self, width: int, height: int):
        """
        Parameters:
            width (int) :
            height (int):
                Dimensions of the aquifer

            >>> g = Grid(8, 8)
            Academic license - for non-commercial use only
        """
        self._width = width
        self._height = height
        self.Model = gu.Model()  # pylint: disable=E1101
        self.wells = _make_zero_tuple_dict(self._width + 1, self._height + 1)
        self.LHS = _make_zero_tuple_dict(self._width + 1, self._height + 1)
        self.well_coords = []
        self.contam_coords = []
        self._bounds = {
            "top": [(i, 0) for i in range(1, self._width + 1)],
            "bottom": [(i, self._height + 1) for i in range(1, self._width + 1)],
            "left": [(0, j) for j in range(1, self._height + 1)],
            "right": [(self._width + 1, j) for j in range(1, self._height + 1)],
        }
        self._heads = self.Model.addVars(
            itertools.product(range(1, self._width + 1), range(1, self._height + 1)),
            name="h",
        )
        self._solved = False
        self._init_boundary_heads()

    def add_boundary_contaminant(self, *args):
        """Adds boundary (left, top, right, bottom) contaminants by
            adding model constraint. Head values within aquifer next to
            boundary are constrained to be greater than the boundary head

            Parameters:
                *args (str): the boundary side ('left', 'top', 'right', 'bottom')

            Example usage:

                >>> g = Grid(8, 8)
                >>> g.add_boundary_contaminant("bottom")
                Traceback (most recent call last):
                    ...
                TypeError: Boundary head must be specified before adding boundary contaminant

                >>> g = Grid(8, 8)
                >>> g.add_boundary_heads(bottom=39) # must be specified before adding the boundary contaminant
                >>> g.add_boundary_contaminant("bottom")
        """
        w, h = self._width, self._height
        for bound in args:
            bound = bound.lower()
            if bound == "right":
                self._boundary_helper(h + 1, w, w + 1, "x")
            elif bound == "left":
                self._boundary_helper(h + 1, 1, 0, "x")
            elif bound == "top":
                self._boundary_helper(w + 1, 1, 0, "y")
            elif bound == "bottom":
                self._boundary_helper(w + 1, h, h + 1, "y")
            else:
                raise ValueError(
                    "Boundary not specified correctly. Should be 'top', 'bottom', 'left', or 'right'"
                )
        self.Model.update()

    def add_boundary_heads(self, **kwargs):
        """Specifies the boundary (left, top, right, bottom) head value.
            If a value is not provided for a boundary, it is assumed to 
            be a no-flow boundary.

            Parameters:
                **kwargs: the boundary side and its head value
                    (ex: Grid.add_boundary_heads(top=29) ).

            Example usage:

                >>> g = Grid(8, 8)
                >>> g.add_boundary_heads(top=29, bottom=39)
        """
        for k in kwargs.keys():
            if k.lower() not in ["top", "bottom", "right", "left"]:
                raise ValueError(
                    "Unknown argument assignment. Must specify head for 'top', 'left', 'right', or 'bottom'."
                )
            self._add_bounds(k.lower(), kwargs[k])

    def add_cell_lengths(self, dx: float, dy: float = None):
        """Specifies the height and width of the cells.

            Parameters:
                dx (float): length in x direction
                dy (float, optional): specifies length in y direction. If not
                    specified, dy = dx.
        """
        self._dx = dx
        self._dy = dy if dy is not None else dx

    def add_contaminant(self, coords: tuple):
        """Adds contaminant to aquifer by adding constrint to model.
            Cell heads surrounding contaminant are greater than or equal
            to contaminant head.

            Parameters:
                coords (tuple): The x,y coordinates of the contaminant.
        """
        self.contam_coords.append(coords)
        x, y = coords
        for i in [-1, 1]:
            if self._heads[x + i, y] is not None:
                self.Model.addConstr(self._heads[x + i, y] >= self._heads[x, y])
            if self._heads[x, y + i] is not None:
                self.Model.addConstr(self._heads[x, y + i] >= self._heads[x, y])
        self.Model.update()

    def add_max_head(self, max_head: float = float("inf")):
        """Adds upper bound for the head values.

            Parameters:
                max_head (float, optional): max head value. Defaults to +inf.
        """
        self._max_head = max_head
        for i in self._heads:
            if isinstance(self._heads[i], gu.Var):  # pylint: disable=E1101
                self._heads[i].setAttr("ub", self._max_head)

    def add_min_head(self, min_head: float = 0):
        """Adds lower bound for the head values.

            Parameters:
                min_head (float, optional): min head value. Defaults to 0.0.
        """
        self._min_head = min_head
        for i in self._heads:
            if isinstance(self._heads[i], gu.Var):  # pylint: disable=E1101
                self._heads[i].setAttr("lb", self._min_head)

    def add_min_pump_rate(self, rate: float):
        """Adds lower bound for the sum of the pumping rates. Adds
            constraint to model.

            Parameters:
                rate (float): min pumping rate.
        """
        self.Model.addConstr(self.wells.sum() >= rate)

    def add_objective(self, num: int):
        """Adds objective function. Select from list:

            1: Maximize total pumping.
            2: Maximize heads throughout aquifer.

            Parameters:
                num (int): The objective function selection.
        """
        if num == 1:
            self.Model.setObjective(
                self.wells.sum(), gu.GRB.MAXIMIZE  # pylint: disable=E1101
            )
        elif num == 2:
            temp = []
            for i in self._heads:
                if isinstance(self._heads[i], gu.Var):  # pylint: disable=E1101
                    temp.append(self._heads[i])
            self.Model.setObjective(sum(temp), gu.GRB.MAXIMIZE)  # pylint: disable=E1101
        else:
            raise ValueError(
                "must add int corresponding to objective number. 1: max pumping. 2: max heads in aquifer."
            )
        self.Model.update()

    def add_trans(self, tx: float, sd: float = 0, distribution: str = "uniform"):
        """Adds transmissivity distribution functions.

            Parameters:
                tx (float): Mean of transmissivity in x-direction.
                sd (float, optional): Standard deviation for distribution 
                    function. Default is 0.0.
                distribution (str): Sets the transmissivity distribution 
                    function. Default is 'uniform.' Select from 'uniform', 
                    'normal' or 'lognormal.'
        """
        self._Tx = self._dist_helper(tx, sd, distribution)
        self._Ty = self._dist_helper(tx, sd, distribution)

    def add_units(self, length: str = "m", time: str = "d"):
        """Adds dimensional units.

            Parameters:
                length (str, optional): units for length. Default to 'm' (meters).
                time (str, optional): units for time. Default to 'd' (days).
        """
        self._length = "m" if length == "m" else "ft"
        self._rate = f"cm{time}" if length == "m" else f"cf{time}"
        self._sink = f"mp{time}" if length == "m" else f"fp{time}"

    def add_wells(self, coords):
        """Adds the wells to the aquifer.

            Parameters:
                coords (List(tuple)): The (x,y) coordinates for the wells. 
                    Must be provided as a single tuple (x,y) or as a list 
                    of tuples [(x1,y1), (x2,y2) ...]
        """
        if type(coords) is tuple:
            coords = [coords]
        if type(coords[0]) is not tuple:
            raise ValueError(
                "Well coordinates must be added as a list of tuples, e.g. [(x1, y1), (x2, y2)]"
            )
        for coord in coords:
            self.wells[coord] = self.Model.addVar(name=f"pump{coord}", lb=float("-inf"))
            self.well_coords.append(coord)
        self.Model.update()

    def optimize(self):
        """This function finalizes the model constraints, updates the model,
            and then tries to solve the model. If optimal solution is found,
            the head values are gathered for printing purposes later.
        """
        self._add_constraints()
        self.Model.update()
        self.Model.optimize()
        if (
            self.Model.Status == 2
        ):  # Status = 2 means optimal solution found and model is solved
            self._solved = True
            self._head_results = self._get_head_results()
            print("OPTIMIZATION COMPLETE")
        else:
            self._head_results = None
            print("OPTIMIZATION NOT COMPLETE")

    def print_grid_results(self, filename: str = "heads"):
        """Prints grid head values to the terminal and to csv file.

            Parameters:
                filename (str, optional): filename (without extension) to 
                    save grid head values.
        """
        if self._solved:
            np.savetxt(f"{filename}.csv", self._head_results, fmt="%.5f", delimiter=",")
            print(f"Results: in {self._length}")
            with np.printoptions(
                precision=3, suppress=True, nanstr="bound", floatmode="fixed"
            ):
                print(self._head_results)
        else:
            print("No results to print.")

    def print_pump_results(self, filename: str = "pumps"):
        """Prints pump values to the terminal and to txt file.

            Parameters:
                filename (str, optional): filename (without extension) to
                    save pump values.
        """
        if self._solved:
            with open(f"{filename}.txt", "w") as f:
                f.write("Pump Results\n")
                print("Pump Results")
                for key in self.wells:
                    if isinstance(self.wells[key], gu.Var):  # pylint: disable=E1101
                        line = f"{key}: {self.wells[key].X:10.5f} {self._sink}, head ({self._length}): {self._heads[key].X:6.2f}"
                        print(line)
                        f.write(line + "\n")
                q = self.wells.sum().getValue()
                print(f"Total Flow (q): {q:.2f} {self._sink}")
                print(f"Total Flow (W): {q*self._dx*self._dy:.2f} {self._rate}")
                f.write(f"Total Flow (q): {q:.2f} {self._sink}\n")
                f.write(f"Total Flow (W): {q*self._dx*self._dy:.2f} {self._rate}\n")
        else:
            print("No results to print.")

    def print_map(self, filename: str = "map", text: bool = False, wells: bool = True):
        """Generates heat map of the aquifer based on head values.

            Parameters:
                filename (str, optional): filename (without extension) to save 
                    image of heat map.
                text (bool, optional): Adds the numerical value to the map.
                    Default is False.
                    Really big aquifers will not show this value very well.
                wells (bool, optional): Adds +/-W to the well locations. 
                    Default is True. The W does not scale with aquifer size, 
                    so it will not show up very well in really large aquifers.
        """
        if self._solved:
            grid = self._head_results
            vmin = np.nanmin(grid)
            vmax = np.nanmax(grid)
            for x, y in self.contam_coords:
                grid[y][x] = vmax + 1
            current_cmap = plt.cm.get_cmap("Blues_r")
            current_cmap.set_bad(color="gray")
            current_cmap.set_over(color="red")
            plt.imshow(grid, cmap=current_cmap, vmin=vmin, vmax=vmax)
            for i, j in self.well_coords:
                label = "-W" if self.wells[(i, j)].X >= 0 else "+W"
                plt.text(i, j, label, ha="center", va="center")
            if text:
                for (j, i), label in np.ndenumerate(grid):
                    plt.text(i, j, round(label, 1), ha="center", va="center")
            plt.colorbar()
            plt.savefig(f"{filename}.pdf", dpi=300)
            plt.show()
        else:
            print("No results to print.")

    def print_report(self):
        """Generates lindo style report."""
        if self._solved:
            lindo(self.Model)
        else:
            print("No results to print.")

    def _add_bounds(self, key, val):
        for coord in self._bounds[key]:
            self._heads[coord] = val

    def _add_constraints(self):
        self._add_lhs()
        self.Model.addConstrs(
            self.LHS[key] == self.wells[key] for key in self.LHS.keys()
        )
        self.Model.update()

    def _add_lhs(self):
        for i in range(1, self._width + 1):
            for j in range(1, self._height + 1):
                tx, ty = next(self._Tx), next(self._Ty)
                self._flow_helper((tx, ty), (i, j), (i, j - 1), "y")  # direction: top
                self._flow_helper((tx, ty), (i, j), (i + 1, j), "x")  # direction: right
                self._flow_helper(
                    (tx, ty), (i, j), (i, j + 1), "y"
                )  # direction: bottom
                self._flow_helper((tx, ty), (i, j), (i - 1, j), "x")  # direction: left
        self.Model.update()

    def _boundary_helper(self, a, b, c, direc):
        if direc == "x":
            if self._heads[c, 1] is None:
                raise TypeError(
                    "Boundary head must be specified before adding boundary contaminant"
                )
            for j in range(1, a):
                self.Model.addConstr(self._heads[b, j] >= self._heads[c, j])
                self.contam_coords.append((c, j))
        elif direc == "y":
            if self._heads[1, c] is None:
                raise TypeError(
                    "Boundary head must be specified before adding boundary contaminant"
                )
            for i in range(1, a):
                self.Model.addConstr(self._heads[i, b] >= self._heads[i, c])
                self.contam_coords.append((i, c))

    @staticmethod
    def _dist_helper(mean, sd, distribution):
        if distribution == "normal" or distribution == "uniform":
            sd = 0 if distribution == "uniform" else sd
            while True:
                yield random.gauss(mean, sd)
        elif "log" in distribution:
            if sd == 0:
                raise ValueError(
                    "Standard deviation for lognormal distribution must be specified."
                )
            while True:
                var = sd / mean
                mu_ln_x = 0.5 * math.log(mean ** 2 / (1 + var ** 2))
                sd_ln_x = math.sqrt(math.log(var ** 2 + 1))
                yield random.lognormvariate(mu_ln_x, sd_ln_x)
        else:
            raise ValueError(
                "Distribution must be specified in add_trans(). Select from 'uniform,' 'normal,' and 'lognormal.'"
            )

    def _flow_helper(self, trans: tuple, coord: tuple, dif: tuple, direc: str):
        tx, ty = trans
        if direc == "x":
            if self._heads[dif] is not None:
                self.LHS[coord] += (
                    tx * (self._heads[dif] - self._heads[coord]) / (self._dx ** 2)
                )
        elif direc == "y":
            if self._heads[dif] is not None:
                self.LHS[coord] += (
                    ty * (self._heads[dif] - self._heads[coord]) / (self._dy ** 2)
                )

    def _get_head_results(self):
        temp = np.full((self._width + 2, self._height + 2), np.nan)
        for i in self._heads:
            if isinstance(self._heads[i], gu.Var):  # pylint: disable=E1101
                temp[i] = self._heads[i].X
            else:
                temp[i] = self._heads[i]
        return temp.T

    def _init_boundary_heads(self):
        for v in self._bounds.values():
            for coord in v:
                self._heads[coord] = None


if __name__ == "__main__":
    import doctest
    doctest.testmod()


    width, height = 8, 8
    g = Grid(width, height)
    g.add_units("ft", "d")
    g.add_min_head(23)
    g.add_cell_lengths(500)
    g.add_trans(12500, sd=500, distribution="normal")

    wells = [
        (3, 3),
        (3, 5),
        (3, 7),
        (4, 3),
        (4, 5),
        (4, 7),
        (5, 3),
        (5, 5),
        (5, 7),
        (6, 3),
        (6, 5),
        (6, 7),
    ]
    g.add_wells(wells)
    g.add_boundary_heads(right=29.4, left=39)
    # g.add_boundary_contaminant("right")
    g.add_contaminant((6, 2))
    g.add_objective(1)

    # g.optimize()
    # g.print_report()
    # g.print_grid_results()
    # g.print_pump_results()
    # g.print_map()
    # print()

"""
from aqpy import aq
g = aq.Grid(8, 8)
g.add_units('ft', 'd')
wells = [(3,3),(5,5),(7,7)]
g.add_wells(wells)
g.add_min_pump_rate(1.5)
g.add_min_head(23)
g.add_cell_lengths(500)
g.add_trans(12500, 100, 'normal')
g.add_boundary_heads(left=39, right=29.4)
g.add_objective(2)
g.optimize()
Optimize a model with 65 rows, 67 columns and 294 nonzeros
Coefficient statistics:
  Matrix range     [5e-02, 1e+00]
  Objective range  [1e+00, 1e+00]
  Bounds range     [2e+01, 2e+01]
  RHS range        [1e+00, 2e+00]
Presolve removed 37 rows and 37 columns
Presolve time: 0.00s
Presolved: 28 rows, 30 columns, 246 nonzeros

Iteration    Objective       Primal Inf.    Dual Inf.      Time
       0    3.7793103e+32   1.140819e+31   3.779310e+02      0s
      29    1.9593423e+03   0.000000e+00   0.000000e+00      0s

Solved in 29 iterations and 0.00 seconds
Optimal objective  1.959342309e+03
OPTIMIZATION IS COMPLETE
>>> g.print_pump_results()
Pump Results
(3, 3):    1.17494 fpd, head (ft):  23.00
(5, 5):   -0.45900 fpd, head (ft):  33.04
(7, 7):    0.78406 fpd, head (ft):  23.00
Total Flow (q): 1.50 fpd
Total Flow (W): 375000.00 cfd
>>> g.print_grid()
Results: in ft
[[ bound  bound  bound  bound  bound  bound  bound  bound  bound  bound]
 [39.000 35.566 32.373 29.939 29.182 29.114 29.194 29.269 29.335 29.400]
 [39.000 35.325 31.615 28.261 28.492 28.968 29.197 29.280 29.335 29.400]
 [39.000 35.121 30.499 23.000 27.557 29.067 29.348 29.318 29.324 29.400]
 [39.000 35.658 32.260 29.181 29.670 30.396 29.808 29.321 29.242 29.400]
 [39.000 36.253 33.700 31.796 31.544 33.040 30.168 28.914 28.925 29.400]
 [39.000 36.653 34.493 32.759 31.672 30.869 28.912 27.243 28.144 29.400]
 [39.000 36.866 34.858 33.076 31.514 29.855 27.367 23.000 27.009 29.400]
 [39.000 36.954 34.995 33.175 31.452 29.669 27.700 26.063 27.491 29.400]
 [ bound  bound  bound  bound  bound  bound  bound  bound  bound  bound]]

"""
