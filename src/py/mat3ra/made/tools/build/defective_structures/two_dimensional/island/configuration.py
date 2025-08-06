from mat3ra.esse.models.materials_category.defective_structures.two_dimensional.island.configuration import (
    IslandDefectConfigurationSchema,
)

from .....build_components.entities.reusable.two_dimensional.slab_stack.configuration import SlabStackConfiguration


class IslandDefectConfiguration(SlabStackConfiguration, IslandDefectConfigurationSchema):
    """
    Configuration for creating an island defect on a slab surface.

    Args:
        stack_components: List containing [slab, isolated_island, vacuum].
    """

    type: str = "IslandDefectConfiguration"
