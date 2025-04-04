from typing import List, Optional

import numpy as np
from mat3ra.utils.mixins import RoundNumericValuesMixin
from pydantic import BaseModel, Field


class Cell(RoundNumericValuesMixin, BaseModel):
    # TODO: figure out how to use ArrayOf3NumberElementsSchema
    vector1: List[float] = Field(default_factory=lambda: [1.0, 0.0, 0.0])
    vector2: List[float] = Field(default_factory=lambda: [0.0, 1.0, 0.0])
    vector3: List[float] = Field(default_factory=lambda: [0.0, 0.0, 1.0])
    __round_precision__ = 6

    @classmethod
    def from_vectors_array(cls, vectors_array: Optional[List[List[float]]] = None) -> "Cell":
        if vectors_array is None:
            vectors_array = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]

        # Ensure vectors are properly converted to lists of floats
        processed_vectors = []
        for vector in vectors_array:
            processed_vector = [float(v) for v in vector]
            processed_vectors.append(processed_vector)

        return cls(vector1=processed_vectors[0], vector2=processed_vectors[1], vector3=processed_vectors[2])

    @property
    def vectors_as_array(self, skip_rounding=False) -> List[List[float]]:
        if skip_rounding:
            return [self.vector1, self.vector2, self.vector3]
        return self.round_array_or_number([self.vector1, self.vector2, self.vector3])

    def to_list(self, skip_rounding=False) -> List[List[float]]:
        _ = self.round_array_or_number
        return [
            self.vector1 if skip_rounding else _(self.vector1),
            self.vector2 if skip_rounding else _(self.vector2),
            self.vector3 if skip_rounding else _(self.vector3),
        ]

    def clone(self) -> "Cell":
        return self.from_vectors_array(self.vectors_as_array)

    def convert_point_to_cartesian(self, point: List[float]) -> List[float]:
        np_vector = np.array(self.vectors_as_array)
        result_list = np.dot(point, np_vector).tolist()
        return self.round_array_or_number(result_list)

    def convert_point_to_crystal(self, point: List[float]) -> List[float]:
        np_vector = np.array(self.vectors_as_array)
        result_list = np.dot(point, np.linalg.inv(np_vector)).tolist()
        return self.round_array_or_number(result_list)

    @property
    def volume(self) -> float:
        volume = np.linalg.det(np.array(self.vectors_as_array))
        return self.round_array_or_number(volume)
