from typing import Type, Union

from mat3ra.made.material import Material
from .configuration import PairDefectConfiguration
from ..point_defect.base.builder import PointDefectBuilder
from ..point_defect.vacancy.builder import VacancyDefectBuilder
from .....build_components import MaterialWithBuildMetadata, TypeConfiguration


class PairDefectBuilder(VacancyDefectBuilder, PointDefectBuilder):
    _ConfigurationType: Type[PairDefectConfiguration] = PairDefectConfiguration

    def _update_material_name(
        self, material: Union[Material, MaterialWithBuildMetadata], configuration: TypeConfiguration
    ) -> Material:
        host_material = configuration.merge_components[0]

        primary_defect_name = (
            configuration.merge_components[1].element.chemical_element.value
            if hasattr(configuration.merge_components[1].element.chemical_element, "value")
            else configuration.merge_components[1].element.chemical_element
        )
        secondary_defect_name = (
            configuration.merge_components[2].element.chemical_element.value
            if hasattr(configuration.merge_components[2].element.chemical_element, "value")
            else configuration.merge_components[2].element.chemical_element
        )
        if host_material:
            material.name = f"{host_material.name}, Pair Defect: {primary_defect_name} + {secondary_defect_name}"

        return material
