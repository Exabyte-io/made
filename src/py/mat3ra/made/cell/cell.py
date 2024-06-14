from typing import List

import numpy as np


class Cell:
    tolerance = 1

    def __init__(self, nested_array):
        self.vector1, self.vector2, self.vector3 = nested_array

    @property
    def vectors_as_array(self) -> List[List[float]]:
        return [[0 if abs(c) < self.tolerance else c for c in v] for v in [self.vector1, self.vector2, self.vector3]]

    def clone(self):
        return self.__class__(self.vectors_as_array)

    def clone_and_scale_by_matrix(self, matrix):
        newCell = self.clone()
        newCell.scale_by_matrix(matrix)
        return newCell

    def convert_point_to_cartesian(self, point):
        np_vector = np.array(self.vectors_as_array)
        return np.dot(point, np_vector)

    def convert_point_to_fractional(self, point):
        np_vector = np.array(self.vectors_as_array)
        return np.dot(point, np.linalg.inv(np_vector))

    def scale_by_matrix(self, matrix):
        np_vector = np.array(self.vectors_as_array)
        self.vector1, self.vector2, self.vector3 = np.dot(matrix, np_vector).tolist()
