from typing import List, Optional

from mat3ra.made.tools.build_components.entities.reusable.base_builder import BaseBuilderParameters


class TerraceBuildParameters(BaseBuilderParameters):
    rotate_to_match_pbc: bool = True
    axis: List[int] = [0, 1, 0]
    angle: float = 0.0
    new_lattice_vectors: Optional[List[List[float]]] = None
