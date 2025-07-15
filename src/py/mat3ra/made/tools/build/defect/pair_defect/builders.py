from typing import Any, Type

from mat3ra.made.material import Material
from mat3ra.made.tools.build import MaterialWithBuildMetadata
from mat3ra.made.tools.build.merge.builders import MergeBuilder
from mat3ra.made.tools.build.defect.point.builders import PointDefectBuilder
from mat3ra.made.tools.build.defect.point.configuration import PointDefectConfiguration
from .configuration import PairDefectConfiguration


class PairDefectBuilder(MergeBuilder):
    """
    Builder for creating pair defects by handling two defect placements.
    Extends MergeBuilder and uses existing coordinate resolution methods.
    """

    _ConfigurationType: Type[PairDefectConfiguration] = PairDefectConfiguration

    def _merge_component_to_material(self, configuration_or_material: Any) -> Material:
        if isinstance(configuration_or_material, PointDefectConfiguration):
            return PointDefectBuilder().get_material(configuration_or_material)
        return super()._merge_component_to_material(configuration_or_material)

    def _update_material_name(self, material: Material, configuration: _ConfigurationType) -> Material:
        host_material = configuration.merge_components[0]

        primary_defect_name = configuration.merge_components[1].type
        secondary_defect_name = configuration.merge_components[2].type
        if host_material:
            material.name = f"{host_material.name}, Pair Defect, {primary_defect_name} + {secondary_defect_name}"

        return material
