import math
from typing import List, Optional

import numpy as np
from mat3ra.code.entity import InMemoryEntityPydantic
from mat3ra.esse.models.properties_directory.structural.lattice.lattice_bravais import (
    LatticeImplicitSchema as LatticeBravaisSchema,
)
from mat3ra.esse.models.properties_directory.structural.lattice.lattice_bravais import (
    LatticeTypeEnum,
    LatticeUnitsSchema,
)
from mat3ra.utils.mixins import RoundNumericValuesMixin

from .cell import Cell

COORDINATE_TOLERANCE = 6


class LatticeVectors(Cell):
    pass


class Lattice(RoundNumericValuesMixin, LatticeBravaisSchema, InMemoryEntityPydantic):
    __types__ = LatticeTypeEnum
    __type_default__ = LatticeBravaisSchema.model_fields["type"].default
    __units_default__ = LatticeBravaisSchema.model_fields["units"].default_factory()

    a: float = 1.0
    b: float = a
    c: float = a
    alpha: float = 90.0
    beta: float = 90.0
    gamma: float = 90.0

    @property
    def vectors(self) -> LatticeVectors:
        vectors = self.calculate_vectors()
        return LatticeVectors.from_vectors_array(vectors)

    def calculate_vectors(self):
        a = self.a
        b = self.b
        c = self.c
        # Convert degrees to radians for trigonometric functions
        alpha_rad = math.radians(self.alpha)
        beta_rad = math.radians(self.beta)
        gamma_rad = math.radians(self.gamma)

        # Calculate cosines and sines of the angles
        cos_alpha = math.cos(alpha_rad)
        cos_beta = math.cos(beta_rad)
        cos_gamma = math.cos(gamma_rad)
        sin_alpha = math.sin(alpha_rad)
        sin_beta = math.sin(beta_rad)

        # Compute gamma star (used in matrix calculation)
        gamma_star = math.acos((cos_alpha * cos_beta - cos_gamma) / (sin_alpha * sin_beta))
        cos_gamma_star = math.cos(gamma_star)
        sin_gamma_star = math.sin(gamma_star)

        # Calculate the vectors
        vector_a = [a * sin_beta, 0.0, a * cos_beta]
        vector_b = [-b * sin_alpha * cos_gamma_star, b * sin_alpha * sin_gamma_star, b * cos_alpha]
        vector_c = [0.0, 0.0, c]

        return [vector_a, vector_b, vector_c]

    @classmethod
    def from_vectors_array(
        cls,
        vectors: List[List[float]],
        units: Optional[LatticeUnitsSchema] = __units_default__,
        type: Optional[LatticeTypeEnum] = __type_default__,
    ) -> "Lattice":
        a = np.linalg.norm(vectors[0])
        b = np.linalg.norm(vectors[1])
        c = np.linalg.norm(vectors[2])
        alpha = np.degrees(np.arccos(np.dot(vectors[1], vectors[2]) / (b * c)))
        beta = np.degrees(np.arccos(np.dot(vectors[0], vectors[2]) / (a * c)))
        gamma = np.degrees(np.arccos(np.dot(vectors[0], vectors[1]) / (a * b)))

        return cls(
            a=float(a),
            b=float(b),
            c=float(c),
            alpha=alpha,
            beta=beta,
            gamma=gamma,
            units=units,
            type=type,
        )

    @property
    def vector_arrays(self) -> List[List[float]]:
        return self.vectors.vector_arrays

    @property
    def vector_arrays_rounded(self) -> List[List[float]]:
        return self.vectors.vector_arrays_rounded

    @property
    def cell_volume(self) -> float:
        return self.vectors.volume

    @property
    def cell_volume_rounded(self) -> float:
        return self.vectors.volume_rounded

    def get_scaled_by_matrix(self, matrix: List[List[float]]):
        """
        Scale the lattice by a matrix.
        Args:
            matrix (List[List[float]]): A 3x3 matrix.
        """
        np_vectors = np.array(self.vector_arrays)
        np_matrix = np.array(matrix)
        scaled_vectors = np.dot(np_matrix, np_vectors).tolist()
        new_lattice = self.from_vectors_array(vectors=scaled_vectors, units=self.units, type=self.type)
        return new_lattice
