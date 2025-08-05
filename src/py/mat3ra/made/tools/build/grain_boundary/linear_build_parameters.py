from typing import Optional

from .build_parameters import GrainBoundaryBuilderParameters


class GrainBoundaryLinearBuilderParameters(GrainBoundaryBuilderParameters):
    max_supercell_matrix_int: Optional[int] = None
    limit_max_int: Optional[int] = 30
    angle_tolerance: float = 0.1
    return_first_match: bool = False
    edge_inclusion_tolerance: float = 1.0
    distance_tolerance: float = 1.0
