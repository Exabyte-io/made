from .configuration import SlabStackConfiguration
from ... import MaterialWithBuildMetadata
from ...slab.builders import SlabBuilder
from ...slab.configurations import SlabConfiguration
from ...stack.builders import StackNComponentsBuilder


class SlabStackBuilder(StackNComponentsBuilder):
    _ConfigurationType = SlabStackConfiguration

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
        self, material: MaterialWithBuildMetadata, configuration: _ConfigurationType
    ) -> MaterialWithBuildMetadata:
        name_suffix = self.get_name_suffix(configuration) or self.name_suffix
        base_name = material.name
        name = f"{base_name}, {name_suffix}"
        material.name = name
        return material
