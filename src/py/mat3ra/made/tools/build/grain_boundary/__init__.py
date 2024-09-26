from typing import Optional

from mat3ra.made.material import Material

from .builders import SlabGrainBoundaryBuilder
from .configuration import SlabGrainBoundaryConfiguration


def create_grain_boundary(
    configuration: SlabGrainBoundaryConfiguration,
    builder: Optional[SlabGrainBoundaryBuilder] = None,
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
