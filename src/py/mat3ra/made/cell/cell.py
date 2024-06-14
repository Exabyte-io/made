import numpy as np


class Cell:
    tolerance = 1

    def __init__(self, nestedArray):
        self.vector1, self.vector2, self.vector3 = nestedArray

    @property
    def vectorsAsArray(self):
        return [
            [0 if abs(c) < self.tolerance else c for c in v]
            for v in [self.vector1, self.vector2, self.vector3]
        ]

    def clone(self):
        return self.__class__(self.vectorsAsArray)

    def cloneAndScaleByMatrix(self, matrix):
        newCell = self.clone()
        newCell.scaleByMatrix(matrix)
        return newCell

    def convertPointToCartesian(self, point):
        return np.dot(point, self.vectorsAsArray)

    def convertPointToFractional(self, point):
        return np.dot(point, np.linalg.inv(self.vectorsAsArray))

    def isPointInsideCell(self, point, tolerance=1):
        fractional_point = self.convertPointToFractional(point)
        return all(0 <= c <= 1 for c in fractional_point)

    def getMostCollinearVectorIndex(self, testVector):
        angles = [np.arccos(np.dot(v, testVector) / (np.linalg.norm(v) * np.linalg.norm(testVector))) for v in self.vectorsAsArray]
        return angles.index(min(angles))

    def scaleByMatrix(self, matrix):
        self.vector1, self.vector2, self.vector3 = np.dot(matrix, self.vectorsAsArray).tolist()


