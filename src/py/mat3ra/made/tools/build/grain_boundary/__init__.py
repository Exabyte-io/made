from typing import Union

from .configuration import SurfaceGrainBoundaryConfiguration, SlabGrainBoundaryConfiguration
from .builders import (
    SurfaceGrainBoundaryBuilder,
    SurfaceGrainBoundaryBuilderParameters,
    SlabGrainBoundaryBuilder,
    SlabGrainBoundaryBuilderParameters,
)


def create_grain_boundary(
    configuration: Union[SurfaceGrainBoundaryConfiguration, SlabGrainBoundaryConfiguration],
    builder_parameters: Union[SurfaceGrainBoundaryBuilderParameters, SlabGrainBoundaryBuilderParameters, None] = None,
):
    """
    Create a grain boundary between two surface phases.

    Args:
        configuration: The configuration of the grain boundary to be created.
        builder_parameters: The parameters to be used by the grain boundary builder.

    Returns:
        The material with the grain boundary added.
    """
    builder = SurfaceGrainBoundaryBuilder(build_parameters=builder_parameters)
    return builder.get_materials(configuration)
