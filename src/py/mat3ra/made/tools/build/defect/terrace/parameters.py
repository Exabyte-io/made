from typing import List

from mat3ra.made.tools.build import BaseBuilderParameters
from mat3ra.made.tools.build.slab.entities import MillerIndices


class TerraceBuildParameters(BaseBuilderParameters):
    rotate_to_match_pbc: bool = True
    axis: List[int] = [0, 1, 0]
    angle: float = 0.0
