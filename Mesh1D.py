"""A simple class that describes a 1-D mesh"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad

from Util import SymbolAsFunc

class Mesh1D:
    """ Mesh1D class definition """
    NodeCoords = []
    ElemConnect = []
    NodeDof = []
    Addressing = []
    LowerLimit = 0
    UpperLimit = 1
    NumElem = 0
    NumDofs = 0
    PolyDegree = 0

    def __init__(self, Boundary):
        self.LowerLimit = min(Boundary)
        self.UpperLimit = max(Boundary)
        self.NodeCoords = [self.LowerLimit, self.UpperLimit]
        self.ElemConnect = [[0, 1]]
        self.NumElem = 1
        self.PolyDegree = 1
        self.InitDofs()

    def print(self):
        """Print info about the mesh"""
        print("Coordinates: ", self.NodeCoords)
        print("Elements: ", self.ElemConnect)

    def SetNumElements(self, NumElements):
        """Divides the domain into elements of equal size"""
        self.NumElem = NumElements
        elemSize = (self.UpperLimit - self.LowerLimit) / NumElements
        self.NodeCoords = [self.LowerLimit + i * elemSize for i in range(0, NumElements + 1)]
        self.ElemConnect = [[i, i+1] for i in range(0, NumElements)]

    def Display(self):
        """Plot the mesh"""
        yaxis = [1 for i in self.NodeCoords]
        plt.plot(self.NodeCoords, yaxis, '-o')
        plt.show()

    def InitDofs(self, polyDegree = 1):
        """ Number the DOFs, creates/erases nodes depending on the polynomial degree """
        self.PolyDegree = polyDegree
        NumNodes = self.NumElem * polyDegree + 1
        self.NumDofs = NumNodes
        intervalSize = (self.UpperLimit - self.LowerLimit) / (NumNodes - 1)
        self.NodeCoords = [self.LowerLimit + i * intervalSize for i in range(0, NumNodes)]
        self.ElemConnect = [[i+j for j in range(0, polyDegree + 1)] for i in range(0, NumNodes-1, polyDegree)]
        self.NodeDof = [i for i in range(0, NumNodes - 2)]
        self.NodeDof.insert(0, NumNodes - 2)
        self.NodeDof.append(NumNodes - 1)
        self.Addressing = [[self.NodeDof[j] for j in self.ElemConnect[i]] for i in range(0, self.NumElem)]

    def GetElemSize(self, elem):
        elc = self.ElemConnect[elem]
        return self.NodeCoords[elc[-1]] - self.NodeCoords[elc[0]]

    def GetAddressingVector(self, elem):
        return self.Addressing[elem]

    def GetElemMatrix(self, elem, k, f):
        """ Computes an elementary matrix """
        # **** Need to do proper integration ****
        m_size = range(self.PolyDegree + 1)
        a_elem = np.matrix([[0.0 for j in m_size] for i in m_size])
        b_elem = np.array([0.0 for i in m_size])
        if self.PolyDegree == 1:
            integrand = SymbolAsFunc(k)
            integralVal = quad(integrand, -1., 1.)[0]
            h = self.GetElemSize(elem)
            a_elem[0,0] = 0.5 / h * integralVal
            a_elem[0,1] = -0.5 / h * integralVal
            a_elem[1,0] = a_elem[0,1]
            a_elem[1,1] = a_elem[0,0]
            b_elem[0] = quad(SymbolAsFunc(f), -1., 1.)[0] * h * 0.5
            b_elem[1] = b_elem[0]

        return a_elem,b_elem

    def ApplyBoundaryCond(self, A, B, k, f, c1, c2):
        A[:,-1] = np.matrix([[0.] for i in range(A.shape[0])])
        A[:,-2] = np.matrix([[0.] for i in range(A.shape[0])])
        A[-1,:] = np.matrix([0. for i in range(A.shape[0])])
        A[-2,:] = np.matrix([0. for i in range(A.shape[0])])
        A[-1,-1] = 1
        A[-2,-2] = 1
        B[-1] = c2
        B[-2] = c1
        # *** Need to get the right index in elementary matrix! ***
        B[self.Addressing[0][1]] = \
            B[self.Addressing[0][1]] - self.GetElemMatrix(0, k, f)[0][1,0] * c1
        B[self.Addressing[-1][-2]] = \
            B[self.Addressing[-1][-2]] - self.GetElemMatrix(self.NumElem - 1, k, f)[0][0,1] * c2


        return A, B
