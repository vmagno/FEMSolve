"""Description of a differential equation solver"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
from sympy import Symbol, diff

from Mesh1D import Mesh1D
from Util import SymbolAsFunc, x

class DifferentialEquation:
    """Differential Equation definition (only heat equation for now)"""

    k = 1
    f = 1
    BCond1 = 1
    BCond2 = 1
    Degree = 5
    ExactSolution = 0

    Solution = None
    LiftFunction = None
    Mesh = None

    def __init__(self, k, f, low, high, bclow, bchigh, exactSolution = 0):
        self.k = k + x - x
        self.f = f + x - x
        self.Mesh = Mesh1D([low, high])
        self.BCond1 = bclow
        self.BCond2 = bchigh

        self.ExactSolution = exactSolution
        if exactSolution == 0:
            self.ExactSolution = symbol - symbol

    def ShowMesh(self):
        self.Mesh.Display()

    def PrintInfo(self):
        print('k = ', self.k)
        print('f = ', self.f)
        print('Boundaries: ', self.Mesh.LowerLimit, ' and ', self.Mesh.UpperLimit)
        print('Bd conditions: ', self.BCond1, ' and ', self.BCond2)

#    def SymbolAsFunc(self, symbol):
#        def _function(a):
#            return symbol.subs(self.x, a)
#        return _function

    def Shape(self, i):
        return (x - self.Mesh.LowerLimit)**i * (x - self.Mesh.UpperLimit)

    def SolveRitz(self, degree):
        #Ignores BC for now
        #shapeFunctions = [self.MakeShapeFunction(i) for i in range(1, Degree + 1)]
        self.Degree = degree
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
            return SymbolAsFunc(diff(shape(i), x) * diff(shape(j), x) * k)
        def integrandB(i):
            return SymbolAsFunc(shape(i) * f - diff(shape(i)) * k * diff(self.LiftFunction))

        A = np.matrix( [[quad(integrandA(i, j), lower, upper)[0] for j in loop]  for i in loop ])
        B = np.array( [quad(integrandB(i), lower, upper)[0] for i in loop] )

        print(A)
        print(B)

        coefs = np.linalg.solve(A, B)
        print(coefs)

        self.Solution = sum([self.Shape(i+1) * coefs[i] for i in range(len(coefs))]) + self.LiftFunction

    def SolveFE(self):
        self.Mesh.SetNumElements(2)
        self.Mesh.InitDofs(1)
        self.PrintInfo()
        print(self.Mesh.NodeCoords)
        print(self.Mesh.ElemConnect)
        print(self.Mesh.NodeDof)
        print(self.Mesh.Addressing)

        m_size = range(self.Mesh.NumDofs)
        A = np.matrix([[0. for j in m_size] for i in m_size])
        B = np.array([0. for i in m_size])

        # Per-element matrix assembly
        for e in range(self.Mesh.NumElem):
            elMat, elB = self.Mesh.GetElemMatrix(e, self.k, self.f)
            addr = self.Mesh.GetAddressingVector(e)
            for i in range(elMat.shape[0]):
                for j in range(elMat.shape[1]):
                    A[addr[i],addr[j]] = A[addr[i],addr[j]] + elMat[i,j]
                    B[addr[i]] = B[addr[i]] + elB[i]
            #print('Element ', i, ': ', elMat)

        print(A)
        print(B)

        # Boundary conditions
        A, B = self.Mesh.ApplyBoundaryCond(A, B, self.k, self.f, self.BCond1, self.BCond2)
        print(A)
        print(B)

        coefs = np.linalg.solve(A, B)
        print(coefs)

        self.Mesh.Display()

    def DisplaySolution(self, degree = 0):
        """ Plots the solution. Solves the equation with the specified degree if
            not done already """
        if degree > 0:
            self.SolveRitz(degree)
        elif self.Solution is None:
            self.SolveRitz(1)

        lower = self.Mesh.LowerLimit
        upper = self.Mesh.UpperLimit

        NumPts = 25
        Interval = (upper - lower) / (NumPts - 1)
        xval = [lower + i * Interval for i in range(NumPts)]
        yval = [SymbolAsFunc(self.Solution)(x) for x in xval]
        yval2 = [SymbolAsFunc(self.LiftFunction)(x) for x in xval]
        yval4 = [SymbolAsFunc(self.ExactSolution)(x) for x in xval]
        plt.plot(xval, yval, '-o', label='Computed solution')
#        plt.plot(xval, yval2, 'r-*', label='Lift function')
        if self.ExactSolution != 0:
            plt.plot(xval, yval4, 'g-^', label='Exact solution')

        plt.legend()
        plt.show()
        plt.gcf().clear()
