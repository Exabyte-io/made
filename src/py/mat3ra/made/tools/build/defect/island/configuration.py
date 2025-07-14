from ..slab.configuration import SlabStackConfiguration


class IslandDefectConfiguration(SlabStackConfiguration):
    """
    Configuration for creating an island defect on a slab surface.

    Args:
        stack_components: List containing [slab, isolated_island, vacuum].
    """

    type: str = "IslandDefectConfiguration"
