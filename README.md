# AqPy

This module is used to solve 2-D aquifer pumping optimization problems.


### Installation

Anaconda is the recommended Python distribution for this module. The Anaconda distribution can be obtained [here](https://www.anaconda.com/distribution/ "Anaconda Distribution Download").


##### Dependencies:

* Gurobi Optimizer (gurobipy)
* matplotlib
* numpy

##### Gurobi Installation:

Once Anaconda has been installed, use the following `conda` commands in your terminal to install Gurobi Optimizer.
```bash
conda config --add channels http://conda.anaconda.org/gurobi

conda install gurobi
```
Gurobi requires a free academic license which can be acquired from [here](https://user.gurobi.com/download/licenses/free-academic "Retrieving a Free Academic license"). Installation and verification of the license is only possible on the university network.

Following the installation of Gurobi, the AqPy package can be installed with the following command:
```bash
pip install git+https://github.com/maskarb/aqpy
```

##### Brief Tutorial:
[ipython notebook](https://github.com/maskarb/AqPy/blob/master/tutorial/tutorial.ipynb)

