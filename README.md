# AqPy

This module is used to solve 2-D aquifer pumping optimization problems.


###Installation
Anaconda is the recommended Python distribution for this module. The Anaconda distribution can be obtained [here](https://www.anaconda.com/distribution/ "Anaconda Distribution Download").


#####Dependencies:
* Gurobi Optimizer (gurobipy)
* matplotlib
* numpy

#####Gurobi Installation:
Once Anaconda has been installed, use the following `conda` commands in your terminal to install Gurobi Optimizer.
```bash
conda config --add channels http://conda.anaconda.org/gurobi

conda install gurobi
```
Gurobi requires a free academic license which can be acquired from [here](https://user.gurobi.com/download/licenses/free-academic "Retrieving a Free Academic license"). Installation and verification of the license is only possible on the university network.




```python
import aqpy
g = aqpy.Grid(8, 8)
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

g.print_pump_results()

Pump Results
(3, 3):    1.17494 fpd, head (ft):  23.00
(5, 5):   -0.45900 fpd, head (ft):  33.04
(7, 7):    0.78406 fpd, head (ft):  23.00
Total Flow (q): 1.50 fpd
Total Flow (W): 375000.00 cfd

g.print_grid()

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
```

