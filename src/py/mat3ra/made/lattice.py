import math
from typing import Any, Dict, List

import numpy as np
from mat3ra.utils.mixins import RoundNumericValuesMixin
from pydantic import BaseModel

from .cell import Cell

HASH_TOLERANCE = 3


class Lattice(RoundNumericValuesMixin, BaseModel):
    a: float = 1.0
    b: float = a
    c: float = a
    alpha: float = 90.0
    beta: float = 90.0
    gamma: float = 90.0
    units: Dict[str, str] = {
        "length": "angstrom",
        "angle": "degree",
    }
    type: str = "TRI"

    @property
    def vectors(self) -> List[List[float]]:
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

        # Return the lattice matrix using the derived trigonometric values
        return [
            [a * sin_beta, 0.0, a * cos_beta],
            [-b * sin_alpha * cos_gamma_star, b * sin_alpha * sin_gamma_star, b * cos_alpha],
            [0.0, 0.0, c],
        ]

    def to_json(self, skip_rounding: bool = False) -> Dict[str, Any]:
        __round__ = RoundNumericValuesMixin.round_array_or_number
        round_func = __round__ if not skip_rounding else lambda x: x
        return {
            "a": round_func(self.a),
            "b": round_func(self.b),
            "c": round_func(self.c),
            "alpha": round_func(self.alpha),
            "beta": round_func(self.beta),
            "gamma": round_func(self.gamma),
            "units": self.units,
            "type": self.type,
            "vectors": self.vectors,
        }

    def clone(self, extra_context: Dict[str, Any]) -> "Lattice":
        return Lattice(**{**self.to_json(), **extra_context})

    @property
    def vector_arrays(self) -> List[List[float]]:
        return self.vectors

    @property
    def cell(self) -> Cell:
        return Cell.from_nested_array(self.vector_arrays)

    def volume(self) -> float:
        np_vector = np.array(self.vector_arrays)
        return abs(np.linalg.det(np_vector))
