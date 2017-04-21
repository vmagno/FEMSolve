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

    SystemMatrix = None
    LoadVector = None
    Coefs = None

    def __init__(self, Boundary):
        self.LowerLimit = min(Boundary)
        self.UpperLimit = max(Boundary)
        self.NodeCoords = [self.LowerLimit, self.UpperLimit]
        self.ElemConnect = [[0, 1]]
        self.NumElem = 1
        self.PolyDegree = 1
        self.InitDofs()

    def PrintInfo(self):
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
            a_elem[0,0] = k / h
            a_elem[0,1] = -k / h
            a_elem[1,0] = a_elem[0,1]
            a_elem[1,1] = a_elem[0,0]
            b_elem[0] = f * h / 4. # ** Why does this work???
            b_elem[1] = b_elem[0]

        return a_elem,b_elem

    def AssembleMatrix(self, k, f):
        m_size = range(self.NumDofs)
        self.SystemMatrix = np.matrix([[0. for j in m_size] for i in m_size])
        self.LoadVector = np.array([0. for i in m_size])
        A = self.SystemMatrix
        B = self.LoadVector

        for e in range(self.NumElem):
            elMat, elB = self.GetElemMatrix(e, k, f)

            print('elMat(', e, '):\n', elMat)
            print('elB(', e, '):\n', elB)
            addr = self.GetAddressingVector(e)
            for i in range(elMat.shape[0]):
                for j in range(elMat.shape[1]):
                    A[addr[i],addr[j]] = A[addr[i],addr[j]] + elMat[i,j]
                    B[addr[i]] = B[addr[i]] + elB[i]

        print(self.SystemMatrix)
        print(self.LoadVector)

        return A, B


    def ApplyBoundaryCond(self, k, f, c1, c2):
        A = self.SystemMatrix
        B = self.LoadVector
        A[:,-1] = np.matrix([[0.] for i in range(A.shape[0])])
        A[:,-2] = np.matrix([[0.] for i in range(A.shape[0])])
        A[-1,:] = np.matrix([0. for i in range(A.shape[0])])
        A[-2,:] = np.matrix([0. for i in range(A.shape[0])])
        A[-1,-1] = 1
        A[-2,-2] = 1
        B[-1] = c2
        B[-2] = c1

        term1 = self.GetElemMatrix(0, k, f)[0][1,0] * c1
        term2 = self.GetElemMatrix(self.NumElem - 1, k, f)[0][0,1] * c2
        print('t1 = ', term1, ', t2 = ', term2)
        # *** Need to get the right index in elementary matrix! ***
        B[self.Addressing[0][1]]   = B[self.Addressing[0][1]] - term1
        B[self.Addressing[-1][-2]] = B[self.Addressing[-1][-2]] - term2

        #B[0] = 4.5

        return A, B

    def SolveSystem(self):
        self.Coefs = np.linalg.solve(self.SystemMatrix, self.LoadVector)
        print('A:\n', self.SystemMatrix)
        print('B:\n', self.LoadVector)
        print('Coefs: ', self.Coefs)
        return self.Coefs

    def GetElemAtPoint(self, pt):

        coord = self.NodeCoords
        connec = self.ElemConnect

        for el in range(self.NumElem):
            low = coord[connec[el][0]]
            high = coord[connec[el][-1]]
            if pt >= low and pt <= high:
                return el, low, high

        return None


    def GetSolValue(self, pt):
        #try:
        coefs = self.Coefs
        el, low, high = self.GetElemAtPoint(pt)
        ptref = (pt - low) * 2 / (high - low) - 1
        value = 0.
        dofs = self.Addressing[el]
        value = value + coefs[dofs[0]] * (1 - ptref) / 2.
        value = value + coefs[dofs[-1]] * (1 + ptref) / 2.

        return value

        #except:
        #    print("Could not get solution at point ", pt, "!!!")

    def DisplaySolution(self, ExactSolution = 0):
        lower = self.LowerLimit
        upper = self.UpperLimit
        NumPts = 25
        Interval = (upper - lower) / (NumPts - 1)
        xval = [lower + i * Interval for i in range(NumPts)]
        yval = [self.GetSolValue(x) for x in xval]

        plt.plot(xval, yval, '-o', label='Computed solution')

        if ExactSolution != 0:
            plt.plot(xval, [ExactSolution(x) for x in xval], 'g-^', label='Exact solution')

        plt.legend()
        plt.show()
