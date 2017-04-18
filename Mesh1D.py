"""A simple class that describes a 1-D mesh"""

import matplotlib.pyplot as plt


class Mesh1D:
    """Mesh1D class definition"""
    NodeCoords = []
    ElemConnect = []
    LowerLimit = 0
    UpperLimit = 1
    NumElem = 0

    def __init__(self, Boundary):
        self.LowerLimit = min(Boundary)
        self.UpperLimit = max(Boundary)
        self.NodeCoords = [self.LowerLimit, self.UpperLimit]
        self.ElemConnect = [[0, 1]]
        self.NumElem = 1

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
