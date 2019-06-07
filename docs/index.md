Module aq
=========

Classes
-------

`Grid(width, height)`
:   The Grid class represents the 2-dimensional aquifer. The class is
    instantiated with width and height and the remaining attributes are
    added with adder methods. The order in which attributes are added
    is important. Check the method documentation if errors are raised.
    
    Attributes:
        width (int): Number of grid cells width
        height (int): Number of grid cells height
    
    Parameters:
        width (int) :
        height (int):
            Dimensions of the aquifer
    
        >>> g = Grid(8, 8)
        Academic license - for non-commercial use only

    ### Methods

    `add_boundary_contaminant(self, *args)`
    :   Adds boundary (left, top, right, bottom) contaminants by
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

    `add_boundary_heads(self, **kwargs)`
    :   Specifies the boundary (left, top, right, bottom) head value.
        If a value is not provided for a boundary, it is assumed to 
        be a no-flow boundary.
        
        Parameters:
            **kwargs (str): the boundary side and its head value.
        
        Example usage:
            >>> g = Grid(8, 8)
            >>> g.add_boundary_heads(top=29, bottom=39)

    `add_cell_lengths(self, dx, dy=None)`
    :   Specifies the height and width of the cells.
        
        Parameters:
            dx (float): length in x direction
            dy (float, optional): specifies length in y direction. If not
                specified, dy = dx.

    `add_contaminant(self, coords)`
    :   Adds contaminant to aquifer by adding constrint to model.
        Cell heads surrounding contaminant are greater than or equal
        to contaminant head.
        
        Parameters:
            coords (tuple): The x,y coordinates of the contaminant.

    `add_max_head(self, max_head=inf)`
    :   Adds upper bound for the head values.
        
        Parameters:
            max_head (float, optional): max head value. Defaults to +inf.

    `add_min_head(self, min_head=0)`
    :   Adds lower bound for the head values.
        
        Parameters:
            min_head (float, optional): min head value. Defaults to 0.0.

    `add_min_pump_rate(self, rate)`
    :   Adds lower bound for the sum of the pumping rates. Adds
        constraint to model.
        
        Parameters:
            rate (float): min pumping rate.

    `add_objective(self, num)`
    :   Adds objective function. Select from list:
        
        1: Maximize total pumping.
        2: Maximize heads throughout aquifer.
        
        Parameters:
            num (int): The objective function selection.

    `add_trans(self, tx, sd=0, distribution='uniform')`
    :   Adds transmissivity distribution functions.
        
        Parameters:
            tx (float): Mean of transmissivity in x-direction.
            sd (float, optional): Standard deviation for distribution 
                function. Default is 0.0.
            distribution (str): Sets the transmissivity distribution 
                function. Default is 'uniform.' Select from 'uniform', 
                'normal' or 'lognormal.'

    `add_units(self, length='m', time='d')`
    :   Adds dimensional units.
        
        Parameters:
            length (str, optional): units for length. Default to 'm' (meters).
            time (str, optional): units for time. Default to 'd' (days).

    `add_wells(self, coords)`
    :   Adds the wells to the aquifer.
        
        Parameters:
            coords (List(tuple)): The (x,y) coordinates for the wells. 
                Must be provided as a single tuple (x,y) or as a list 
                of tuples [(x1,y1), (x2,y2) ...]

    `optimize(self)`
    :   This function finalizes the model constraints, updates the model,
        and then tries to solve the model. If optimal solution is found,
        the head values are gathered for printing purposes later.

    `print_grid_results(self, filename='heads')`
    :   Prints grid head values to the terminal and to csv file.
        
        Parameters:
            filename (str, optional): filename (without extension) to 
                save grid head values.

    `print_map(self, filename='map', text=False, wells=True)`
    :   Generates heat map of the aquifer based on head values.
        
        Parameters:
            filename (str, optional): filename (without extension) to save 
                image of heat map.
            text (bool, optional): Adds the numerical value to the map.
                Default is False.
                Really big aquifers will not show this value very well.
            wells (bool, optional): Adds +/-W to the well locations. 
                Default is True. The W does not scale with aquifer size, 
                so it will not show up very well in really large aquifers.

    `print_pump_results(self, filename='pumps')`
    :   Prints pump values to the terminal and to txt file.
        
        Parameters:
            filename (str, optional): filename (without extension) to
                save pump values.

    `print_report(self)`
    :   Generates lindo style report.