from typing import List, Optional

import numpy as np
from mat3ra.utils.mixins import RoundNumericValuesMixin
from pydantic import BaseModel, Field

from mat3ra.code.vector import RoundedVector3D

DEFAULT_CELL = np.eye(3).tolist()

from mat3ra.esse.models.properties_directory.structural.lattice.lattice_vectors import (
    LatticeExplicitUnit as CellSchema,
)


class Cell(RoundNumericValuesMixin, CellSchema):
    a: RoundedVector3D = Field(default_factory=lambda: RoundedVector3D(DEFAULT_CELL[0]))
    b: RoundedVector3D = Field(default_factory=lambda: RoundedVector3D(DEFAULT_CELL[1]))
    c: RoundedVector3D = Field(default_factory=lambda: RoundedVector3D(DEFAULT_CELL[2]))
    __round_precision__ = 6
    __default_vectors__ = DEFAULT_CELL

    @classmethod
    def from_vectors_array(cls, vectors: Optional[List[List[float]]] = None) -> "Cell":
        return cls(
            a=RoundedVector3D(vectors[0]),
            b=RoundedVector3D(vectors[1]),
            c=RoundedVector3D(vectors[2]),
        )

    @property
    def vector_arrays(self, skip_rounding=False) -> List[List[float]]:
        if skip_rounding:
            return [self.a.value, self.b.value, self.c.value]
        return [self.a.value_rounded, self.b.value_rounded, self.c.value_rounded]

    def convert_point_to_cartesian(self, point: List[float]) -> List[float]:
        np_vector = np.array(self.vector_arrays)
        result_list = np.dot(point, np_vector).tolist()
        return self.round_array_or_number(result_list)

    def convert_point_to_crystal(self, point: List[float]) -> List[float]:
        np_vector = np.array(self.vector_arrays)
        result_list = np.dot(point, np.linalg.inv(np_vector)).tolist()
        return self.round_array_or_number(result_list)

    @property
    def volume(self) -> float:
        volume = np.linalg.det(np.array(self.vector_arrays))
        return self.round_array_or_number(volume)
