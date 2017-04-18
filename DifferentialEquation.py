"""Description of a differential equation solver"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
from sympy import Symbol, diff

from Mesh1D import Mesh1D

class DifferentialEquation:
    """Differential Equation definition (only heat equation for now)"""

    x = Symbol('x')
    k = 1+x
    f = x**3 + x
    BCond1 = 1
    BCond2 = 1
    Degree = 5
    ExactSolution = 0

    Solution = None
    LiftFunction = None
    Mesh = None

    def __init__(self, symbol, k, f, low, high, bclow, bchigh, exactSolution = 0):
        self.x = symbol
        self.k = k
        self.f = f
        self.Mesh = Mesh1D([low, high])
        self.BCond1 = bclow
        self.BCond2 = bchigh

        self.ExactSolution = exactSolution
        if exactSolution == 0:
            self.ExactSolution = symbol - symbol

    def ShowMesh(self):
        self.Mesh.Display()

    def SymbolAsFunc(self, symbol):
        def _function(a):
            return symbol.subs(self.x, a)
        return _function

    def Shape(self, i):
        return (self.x - self.Mesh.LowerLimit)**i * (self.x - self.Mesh.UpperLimit)

    def Solve(self, degree):
        #Ignores BC for now
        #shapeFunctions = [self.MakeShapeFunction(i) for i in range(1, Degree + 1)]
        self.Degree = degree
        x = self.x
        shape = self.Shape
        k = self.k
        f = self.f
        lower = self.Mesh.LowerLimit
        upper = self.Mesh.UpperLimit
        loop = range(1, self.Degree + 1)

        slope = (self.BCond2 - self.BCond1) / (upper - lower)
        zero = self.BCond1 - slope * lower
        self.LiftFunction = slope * x + zero

        def integrandA(i, j):
            return self.SymbolAsFunc(diff(shape(i), x) * diff(shape(j), x) * k)
        def integrandB(i):
            return self.SymbolAsFunc(shape(i) * f - diff(shape(i)) * k * diff(self.LiftFunction))

        A = np.matrix( [[quad(integrandA(i, j), lower, upper)[0] for j in loop]  for i in loop ])
        B = np.array( [quad(integrandB(i), lower, upper)[0] for i in loop] )

        print(A)
        print(B)

        coefs = np.linalg.solve(A, B)
        print(coefs)

        self.Solution = sum([self.Shape(i+1) * coefs[i] for i in range(len(coefs))]) + self.LiftFunction

    def DisplaySolution(self, degree = 0):

        if degree > 0:
            self.Solve(degree)
        elif self.Solution is None:
            self.Solve(1)

        lower = self.Mesh.LowerLimit
        upper = self.Mesh.UpperLimit
        T0 = self.BCond1
        T1 = self.BCond2

        k0 = (T1 - T0 + (upper**2 - lower**2) / 2) / (upper - lower)
        k1 = (T0 + lower**2 / 2) - (lower * k0)

        k0_2 = (T1 - T0 + (upper**3 - lower**3) / 6) / (upper - lower)
        k1_2 = (T0 + lower**3 / 6) - (lower * k0_2)

        NumPts = 25
        Interval = (upper - lower) / (NumPts - 1)
        xval = [lower + i * Interval for i in range(NumPts)]
        yval = [self.SymbolAsFunc(self.Solution)(x) for x in xval]
        yval2 = [self.SymbolAsFunc(self.LiftFunction)(x) for x in xval]
        yval4 = [self.SymbolAsFunc(self.ExactSolution)(x) for x in xval]
        plt.plot(xval, yval, '-o', label='Computed solution')
#        plt.plot(xval, yval2, 'r-*', label='Lift function')
        if self.ExactSolution != 0:
            plt.plot(xval, yval4, 'g-^', label='Exact solution')

        plt.legend()
        plt.show()
        plt.gcf().clear()
