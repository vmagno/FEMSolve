"""A simple class that describes a 1-D mesh"""

import matplotlib.pyplot as plt
import numpy as np

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

    def GetElemMatrix(self, elem, k):
        """ Computes an elementary matrix """
        m_size = self.PolyDegree + 1
        a_elem = np.matrix([[0 for j in m_size] for i in m_size])
        b_elem = np.array([0 for i in m_size])
        if self.PolyDegree == 1:
#            a_elem[0][0] =


