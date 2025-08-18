from .configuration import IslandDefectConfiguration
from .....build_components.entities.reusable.two_dimensional.slab_stack.builder import SlabStackBuilder


class IslandDefectBuilder(SlabStackBuilder):
    """
    Builder for creating island defects by merging a slab with a void site.
    The void site defines which atoms to remove from the slab.
    """

    _ConfigurationType = IslandDefectConfiguration

    @property
    def name_suffix(self) -> str:
        return "Island"
