from typing import Dict, Type

from .....build_components import MaterialWithBuildMetadata
from .....build_components.entities.core.zero_dimensional.atom.builder import AtomAtCoordinateBuilder
from .....build_components.entities.core.zero_dimensional.atom.configuration import AtomAtCoordinateConfiguration
from .....build_components.operations.core.combinations.merge import MergeBuilder

from .configuration import PassivationConfiguration


class PassivationBuilder(MergeBuilder):
    @property
    def merge_component_types_conversion_map(self) -> Dict[Type, Type]:
        return {
            AtomAtCoordinateConfiguration: AtomAtCoordinateBuilder,
        }

    def _update_material_name(
        self, material: MaterialWithBuildMetadata, configuration: PassivationConfiguration
    ) -> MaterialWithBuildMetadata:
        material_name = configuration.material.name
        material.name = f"{material_name}, {configuration.passivant}-passivated"
        return material
