from typing import List

from mat3ra.made.tools.build_components.entities.reusable.base_builder import BaseBuilderParameters


class SlabBuilderParameters(BaseBuilderParameters):
    use_orthogonal_c: bool = True
    xy_supercell_matrix: List[List[int]] = [[1, 0], [0, 1]]
