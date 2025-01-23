import json
import math
from typing import Any, Dict, List, Optional

import numpy as np
from mat3ra.utils.mixins import RoundNumericValuesMixin
from pydantic import BaseModel

from .cell import Cell

HASH_TOLERANCE = 3
DEFAULT_UNITS = {"length": "angstrom", "angle": "degree"}
DEFAULT_TYPE = "TRI"


class LatticeVectors(RoundNumericValuesMixin, BaseModel):
    """
    A class to represent the lattice vectors.
    """

    a: List[float] = [1.0, 0.0, 0.0]
    b: List[float] = [0.0, 1.0, 0.0]
    c: List[float] = [0.0, 0.0, 1.0]

    def to_json(self, skip_rounding=False) -> Dict[str, List[float]]:
        json_value = {
            "a": self.round_array_or_number(self.a, skip_rounding) if not skip_rounding else self.a,
            "b": self.round_array_or_number(self.b, skip_rounding) if not skip_rounding else self.b,
            "c": self.round_array_or_number(self.c, skip_rounding) if not skip_rounding else self.c,
        }
        return json.loads(json.dumps(json_value))


class Lattice(RoundNumericValuesMixin, BaseModel):
    a: float = 1.0
    b: float = a
    c: float = a
    alpha: float = 90.0
    beta: float = 90.0
    gamma: float = 90.0
    units: Dict[str, str] = DEFAULT_UNITS
    type: str = DEFAULT_TYPE

    @property
    def vectors(self) -> LatticeVectors:
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

        return LatticeVectors(a=vector_a, b=vector_b, c=vector_c)

    @classmethod
    def from_vectors_array(
        cls, vectors: List[List[float]], units: Optional[Dict[str, str]] = None, type: Optional[str] = None
    ) -> "Lattice":
        """
        Create a Lattice object from a nested array of vectors.
        Args:
            vectors (List[List[float]]): A nested array of vectors.
        Returns:
            Lattice: A Lattice object.
        """
        a = np.linalg.norm(vectors[0])
        b = np.linalg.norm(vectors[1])
        c = np.linalg.norm(vectors[2])
        alpha = np.degrees(np.arccos(np.dot(vectors[1], vectors[2]) / (b * c)))
        beta = np.degrees(np.arccos(np.dot(vectors[0], vectors[2]) / (a * c)))
        gamma = np.degrees(np.arccos(np.dot(vectors[0], vectors[1]) / (a * b)))
        if units is None:
            units = DEFAULT_UNITS
        if type is None:
            type = DEFAULT_TYPE
        return cls(a=float(a), b=float(b), c=float(c), alpha=alpha, beta=beta, gamma=gamma, units=units, type=type)

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
            "vectors": self.vectors.to_json(skip_rounding),
        }

    def clone(self, extra_context: Optional[Dict[str, Any]] = None) -> "Lattice":
        if extra_context is None:
            extra_context = {}
        return Lattice(**{**self.to_json(), **extra_context})

    @property
    def vector_arrays(self) -> List[List[float]]:
        """Returns lattice vectors as a nested array"""
        v = self.vectors
        return [v.a, v.b, v.c]

    @property
    def cell(self) -> Cell:
        return Cell.from_vectors_array(self.vector_arrays)

    def volume(self) -> float:
        np_vector = np.array(self.vector_arrays)
        return abs(np.linalg.det(np_vector))
