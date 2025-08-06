from typing import List

from ...three_dimensional.crystal_lattice_base.build_parameters import (
    BaseBuilderParameters,
)


class SlabBuilderParameters(BaseBuilderParameters):
    use_orthogonal_c: bool = True
    xy_supercell_matrix: List[List[int]] = [[1, 0], [0, 1]]
