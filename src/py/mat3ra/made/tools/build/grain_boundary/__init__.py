from typing import Optional

from mat3ra.made.material import Material

from .builders import GrainBoundaryBuilder
from .configuration import GrainBoundaryConfiguration


def create_grain_boundary(
    configuration: GrainBoundaryConfiguration,
    builder: Optional[GrainBoundaryBuilder] = None,
) -> Material:
    if builder is None:
        builder = GrainBoundaryBuilder()
    return builder.get_material(configuration)
