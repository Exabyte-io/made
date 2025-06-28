from typing import Optional, Union

from mat3ra.made.material import Material

from .builders import GrainBoundaryBuilder, GrainBoundaryWithVacuumBuilder
from .configuration import GrainBoundaryConfiguration, GrainBoundaryWithVacuumConfiguration
from .helpers import (
    create_grain_boundary,
    create_grain_boundary_planar,
    create_grain_boundary_planar_with_vacuum,
)

__all__ = [
    "GrainBoundaryBuilder",
    "GrainBoundaryWithVacuumBuilder",
    "GrainBoundaryConfiguration",
    "GrainBoundaryWithVacuumConfiguration",
    "create_grain_boundary",
    "create_grain_boundary_planar",
    "create_grain_boundary_planar_with_vacuum",
]
