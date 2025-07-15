from typing import Type

from mat3ra.made.material import Material
from mat3ra.made.tools.build.defect.point.builders import (
    PointDefectBuilder,
    VacancyDefectBuilder,
)
from .configuration import PairDefectConfiguration


class PairDefectBuilder(VacancyDefectBuilder, PointDefectBuilder):
    """
    Builder for creating pair defects by handling two defect placements.
    Extends MergeBuilder and uses component type conversion maps.
    """

    _ConfigurationType: Type[PairDefectConfiguration] = PairDefectConfiguration

    def _update_material_name(self, material: Material, configuration: _ConfigurationType) -> Material:
        host_material = configuration.merge_components[0]

        primary_defect_name = configuration.merge_components[1].element.chemical_element.value
        secondary_defect_name = configuration.merge_components[2].element.chemical_element.value
        if host_material:
            material.name = f"{host_material.name}, Pair Defect: {primary_defect_name} + {secondary_defect_name}"

        return material
