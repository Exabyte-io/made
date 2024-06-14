from typing import List, Dict, Any

import math
import numpy as np
from scipy.spatial import ConvexHull

from .cell import Cell
from .constants import HASH_TOLERANCE
from .math import round_value
from .unit_cell import UnitCell
from .utils import get_bravais_lattice_type, get_bravais_lattice_type_extended


class Lattice:
    def __init__(self, config: Dict[str, Any] = {}) -> None:
        self.a = config.get("a", 0)
        self.b = config.get("b", 0)
        self.c = config.get("c", 0)
        self.alpha = config.get("alpha", 0)
        self.beta = config.get("beta", 0)
        self.gamma = config.get("gamma", 0)
        self.units = config.get("units", {"length": "angstrom", "angle": "degree"})
        self.type = config.get("type", "")
        self.vectors = config.get("vectors", [])

    @staticmethod
    def from_vectors(config: Dict[str, Any]) -> "Lattice":
        return Lattice(config)

    def to_json(self, skip_rounding: bool = False) -> Dict[str, Any]:
        round_func = round_value if not skip_rounding else lambda x: x
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
        return Lattice({**self.to_json(), **extra_context})

    @property
    def vector_arrays(self) -> List[List[float]]:
        return self.vectors

    @property
    def cell(self) -> Cell:
        return Cell(self.vector_arrays)

    @property
    def type_label(self) -> str:
        return get_bravais_lattice_type(self.type)

    @property
    def type_extended(self) -> str:
        return get_bravais_lattice_type_extended(
            self.a, self.b, self.c, self.alpha, self.beta, self.gamma, self.type
        )

    @property
    def volume(self) -> float:
        return abs(np.linalg.det(self.vector_arrays))

    @staticmethod
    def get_default_primitive_lattice_config_by_type(
        lattice_config: Dict[str, Any]
    ) -> Dict[str, Any]:
        f_ = round_value
        p_cell = primitive_cell(lattice_config, True)
        new_lattice = Lattice.from_vector_arrays(p_cell, lattice_config["type"])
        k = lattice_config["a"] / new_lattice.a
        return {
            **new_lattice.to_json(),
            "a": f_(new_lattice.a * k),
            "b": f_(new_lattice.b * k),
            "c": f_(new_lattice.c * k),
            "alpha": f_(new_lattice.alpha),
            "beta": f_(new_lattice.beta),
            "gamma": f_(new_lattice.gamma),
        }

    @property
    def unit_cell(self) -> UnitCell:
        vectors = np.array(self.vector_arrays).flatten().tolist() + [self.units["length"]]
        return UnitCell(vectors)

    def get_hash_string(self, is_scaled: bool = False) -> str:
        scale_k = self.a if is_scaled else 1
        scaled_lattice = {
            "a": self.a / scale_k,
            "b": self.b / scale_k,
            "c": self.c / scale_k,
            "alpha": self.alpha,
            "beta": self.beta,
            "gamma": self.gamma,
        }
        return ";".join(
            [
                str(round(x, HASH_TOLERANCE))
                for x in [
                    scaled_lattice["a"],
                    scaled_lattice["b"],
                    scaled_lattice["c"],
                    scaled_lattice["alpha"],
                    scaled_lattice["beta"],
                    scaled_lattice["gamma"],
                ]
            ]
        )


