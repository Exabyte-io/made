from typing import Type, Union

from mat3ra.made.material import Material
from mat3ra.made.tools.build.defect.point.builders import (
    PointDefectBuilder,
    VacancyDefectBuilder,
)
from .configuration import PairDefectConfiguration
from ... import MaterialWithBuildMetadata


class PairDefectBuilder(VacancyDefectBuilder, PointDefectBuilder):
    _ConfigurationType: Type[PairDefectConfiguration] = PairDefectConfiguration

    def _update_material_name(
        self, material: Union[Material, MaterialWithBuildMetadata], configuration: _ConfigurationType
    ) -> Material:
        host_material = configuration.merge_components[0]

        primary_defect_name = configuration.merge_components[1].element.chemical_element.value
        secondary_defect_name = configuration.merge_components[2].element.chemical_element.value
        if host_material:
            material.name = f"{host_material.name}, Pair Defect: {primary_defect_name} + {secondary_defect_name}"

        return material
