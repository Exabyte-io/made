from typing import List

import numpy as np
from mat3ra.esse.models.core.primitive.array_of_3_numbers import ArrayOf3NumberElementsSchema
from mat3ra.utils.mixins import RoundNumericValuesMixin
from pydantic import BaseModel


class Cell(RoundNumericValuesMixin, BaseModel):
    # TODO: figure out how to use
    vector1: ArrayOf3NumberElementsSchema = [1, 0, 0]
    vector2: ArrayOf3NumberElementsSchema = [0, 1, 0]
    vector3: ArrayOf3NumberElementsSchema = [0, 0, 1]
    __round_precision__ = 1e-6

    @classmethod
    def from_nested_array(cls, nested_array):
        if not nested_array:
            nested_array = [cls.vector1, cls.vector2, cls.vector3]
        return cls(vector1=nested_array[0], vector2=nested_array[1], vector3=nested_array[2])

    def __init__(self, vector1=[1, 0, 0], vector2=[0, 1, 0], vector3=[0, 0, 1]):
        super().__init__(**{"vector1": vector1, "vector2": vector2, "vector3": vector3})

    @property
    def vectors_as_array(self, skip_rounding=False) -> List[ArrayOf3NumberElementsSchema]:
        if skip_rounding:
            return [self.vector1, self.vector2, self.vector3]
        return self.round_array_or_number([self.vector1, self.vector2, self.vector3])

    def to_json(self, skip_rounding=False):
        _ = self.round_array_or_number
        if skip_rounding:
            return {
                "vector1": _(self.vector1) if skip_rounding else self.vector1,
                "vector2": _(self.vector2) if skip_rounding else self.vector2,
                "vector3": _(self.vector3) if skip_rounding else self.vector3,
            }

    def clone(self):
        return self.from_nested_array(self.vectors_as_array)

    def clone_and_scale_by_matrix(self, matrix):
        new_cell = self.clone()
        new_cell.scale_by_matrix(matrix)
        return new_cell

    def convert_point_to_cartesian(self, point):
        np_vector = np.array(self.vectors_as_array)
        return np.dot(point, np_vector)

    def convert_point_to_fractional(self, point):
        np_vector = np.array(self.vectors_as_array)
        return np.dot(point, np.linalg.inv(np_vector))

    def scale_by_matrix(self, matrix):
        np_vector = np.array(self.vectors_as_array)
        self.vector1, self.vector2, self.vector3 = np.dot(matrix, np_vector).tolist()
