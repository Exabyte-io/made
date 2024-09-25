from typing import Optional

from mat3ra.made.material import Material

from .builders import GrainBoundaryBuilder
from .configuration import GrainBoundaryConfiguration


def create_grain_boundary(
    configuration: GrainBoundaryConfiguration,
    builder: Optional[GrainBoundaryBuilder] = None,
) -> Material:
    """
    Create a grain boundary according to provided configuration with selected builder.
    Args:
        configuration (GrainBoundaryConfiguration): The configuration of the grain boundary.
        builder (Optional[GrainBoundaryBuilder]): The builder to use for creating the grain boundary.
    Returns:
        Material: The material with the grain boundary.

    """
    if builder is None:
        builder = GrainBoundaryBuilder()
    return builder.get_material(configuration)
