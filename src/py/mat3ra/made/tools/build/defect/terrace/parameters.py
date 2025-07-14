from typing import List

from mat3ra.made.tools.build import BaseBuilderParameters


class TerraceBuildParameters(BaseBuilderParameters):
    rotate_to_match_pbc: bool = True
    axis: List[int] = [0, 1, 0]
    angle: float = 0.0
    new_lattice_vectors: List[List[float]] = None
