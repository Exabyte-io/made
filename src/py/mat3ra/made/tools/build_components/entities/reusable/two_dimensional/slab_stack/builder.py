from typing import Type

from .....build_components import MaterialWithBuildMetadata

from ...slab.slab.builder import SlabBuilder
from ...slab.slab.configuration import SlabConfiguration
from ...stack.builder import StackNComponentsBuilder
from .configuration import SlabStackConfiguration


class SlabStackBuilder(StackNComponentsBuilder):
    _ConfigurationType: Type[SlabStackConfiguration] = SlabStackConfiguration

    @property
    def stack_component_types_conversion_map(self):
        return {
            **super().stack_component_types_conversion_map,
            SlabConfiguration: SlabBuilder,
        }

    @property
    def name_suffix(self) -> str:
        return "Slab Stack"

    def get_name_suffix(self, configuration: SlabStackConfiguration) -> str:
        return self.name_suffix

    def _update_material_name(
        self, material: MaterialWithBuildMetadata, configuration: TypeConfiguration
    ) -> MaterialWithBuildMetadata:
        name_suffix = self.get_name_suffix(configuration) or self.name_suffix
        base_name = material.name
        name = f"{base_name}, {name_suffix}"
        material.name = name
        return material
