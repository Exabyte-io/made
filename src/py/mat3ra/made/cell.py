from typing import List, Optional

import numpy as np
from mat3ra.utils.mixins import RoundNumericValuesMixin
from pydantic import BaseModel, Field


class Cell(RoundNumericValuesMixin, BaseModel):
    # TODO: figure out how to use ArrayOf3NumberElementsSchema
    vector1: List[float] = Field(default_factory=lambda: [1, 0, 0])
    vector2: List[float] = Field(default_factory=lambda: [0, 1, 0])
    vector3: List[float] = Field(default_factory=lambda: [0, 0, 1])
    __round_precision__ = 6

    @classmethod
    def from_vectors_array(cls, vectors_array: Optional[List[List[float]]] = None) -> "Cell":
        if vectors_array is None:
            vectors_array = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        return cls(vector1=vectors_array[0], vector2=vectors_array[1], vector3=vectors_array[2])

    @property
    def vectors_as_array(self, skip_rounding=False) -> List[List[float]]:
        if skip_rounding:
            return [self.vector1, self.vector2, self.vector3]
        return self.round_array_or_number([self.vector1, self.vector2, self.vector3])

    def to_json(self, skip_rounding=False):
        _ = self.round_array_or_number
        return [
            self.vector1 if skip_rounding else _(self.vector1),
            self.vector2 if skip_rounding else _(self.vector2),
            self.vector3 if skip_rounding else _(self.vector3),
        ]

    def clone(self) -> "Cell":
        return self.from_vectors_array(self.vectors_as_array)

    def clone_and_scale_by_matrix(self, matrix: List[List[float]]) -> "Cell":
        new_cell = self.clone()
        new_cell.scale_by_matrix(matrix)
        return new_cell

    def convert_point_to_cartesian(self, point: List[float]) -> List[float]:
        np_vector = np.array(self.vectors_as_array)
        return np.dot(point, np_vector)

    def convert_point_to_crystal(self, point: List[float]) -> List[float]:
        np_vector = np.array(self.vectors_as_array)
        return np.dot(point, np.linalg.inv(np_vector))

    def scale_by_matrix(self, matrix: List[List[float]]):
        np_vector = np.array(self.vectors_as_array)
        self.vector1, self.vector2, self.vector3 = np.dot(np.array(matrix), np_vector).tolist()

    @property
    def volume(self) -> float:
        return np.linalg.det(np.array(self.vectors_as_array))
