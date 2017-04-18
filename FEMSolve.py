#!/usr/bin/python

"""Main FEM solver module."""

from sympy import Symbol

from DifferentialEquation import DifferentialEquation
from Util import x

def MeshTest():
    """Basic function for testing the Mesh!D module"""

    #abc = DifferentialEquation(x, 1, 4*x, 0, 1, 1, 3, -2*x**3 / 3 + 8*x / 3 + 1)
    abc = DifferentialEquation(1, 3*x**2 - x + 1, -1, 2, 2, 0,
                               -(x**4)/4 + (x**3)/6 - (x**2)/2 + 7.0*x/12 + 7.0/2)

    abc.SolveFE()
    abc.DisplaySolution(2)
    #abc.ShowMesh()


if __name__ == "__main__":
    print("Yayy")
    MeshTest()
