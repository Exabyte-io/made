from typing import Dict, Type

from .configuration import PassivationConfiguration
from .. import MaterialWithBuildMetadata
from ..defect.point.builders import AtomAtCoordinateConfiguration, AtomAtCoordinateBuilder
from ..merge import MergeBuilder


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
