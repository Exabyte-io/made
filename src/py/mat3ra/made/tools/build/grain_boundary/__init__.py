from typing import Optional, Union

from mat3ra.made.material import Material

from .builders import SlabGrainBoundaryBuilder, SurfaceGrainBoundaryBuilder, SurfaceGrainBoundaryBuilderParameters
from .configuration import (
    SlabGrainBoundaryConfiguration,
    SurfaceGrainBoundaryConfiguration,
)


def create_grain_boundary(
    configuration: Union[SlabGrainBoundaryConfiguration, SurfaceGrainBoundaryConfiguration],
    builder: Union[SlabGrainBoundaryBuilder, SurfaceGrainBoundaryBuilder, None] = None,
) -> Material:
    """
    Create a grain boundary according to provided configuration with selected builder.
    Args:
        configuration (SlabGrainBoundaryConfiguration): The configuration of the grain boundary.
        builder (Optional[SlabGrainBoundaryBuilder]): The builder to use for creating the grain boundary.
    Returns:
        Material: The material with the grain boundary.

    """
    if builder is None:
        builder = SlabGrainBoundaryBuilder()
    return builder.get_material(configuration)
