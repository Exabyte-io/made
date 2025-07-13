from ..slab.configuration import IslandDefectConfiguration
from ...slab.builders import SlabBuilder
from ...slab.configurations import SlabConfiguration
from ...stack.builders import StackNComponentsBuilder


class SlabStackBuilder(StackNComponentsBuilder):
    @property
    def stack_component_types_conversion_map(self):
        return {
            **super().stack_component_types_conversion_map,
            SlabConfiguration: SlabBuilder,
        }


class IslandDefectBuilder(SlabStackBuilder):
    """
    Builder for creating island defects by merging a slab with a void site.
    The void site defines which atoms to remove from the slab.
    """

    _ConfigurationType = IslandDefectConfiguration
