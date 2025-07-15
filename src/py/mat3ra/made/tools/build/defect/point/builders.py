from typing import Any, Type

from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from mat3ra.esse.models.materials_category_components.entities.core.zero_dimensional.atom import (
    AtomSchema,
)

from mat3ra.made.material import Material
from .configuration import (
    PointDefectConfiguration,
    PointDefectSiteConfiguration,
    VacancyDefectConfiguration,
    SubstitutionalDefectConfiguration,
    InterstitialDefectConfiguration,
)
from ... import BaseSingleBuilder, MaterialWithBuildMetadata
from ...merge.builders import MergeBuilder
from ...vacuum.builders import VacuumBuilder
from ...vacuum.configuration import VacuumConfiguration


class AtomAtCoordinateConfiguration(VacuumConfiguration, PointDefectSiteConfiguration):
    element: AtomSchema


class AtomAtCoordinateBuilder(VacuumBuilder):
    _ConfigurationType = AtomAtCoordinateConfiguration

    def _generate(self, configuration: AtomAtCoordinateConfiguration) -> MaterialWithBuildMetadata:
        vacuum_configuration = VacuumConfiguration(
            crystal=configuration.crystal,
            size=configuration.crystal.lattice.c,
            direction=AxisEnum.z,
        )
        new_material = super()._generate(vacuum_configuration)
        new_material.basis.add_atom(
            element=configuration.element.chemical_element.value,
            coordinate=configuration.coordinate,
        )
        return new_material


class PointDefectSiteBuilder(BaseSingleBuilder):
    """
    Builder class for creating a material from a PointDefectSite configuration.
    """

    _ConfigurationType = PointDefectSiteConfiguration

    def _generate(self, configuration: PointDefectSiteConfiguration) -> MaterialWithBuildMetadata:
        new_material = MaterialWithBuildMetadata.create(
            {
                "name": configuration.crystal.name,
                "lattice": configuration.crystal.lattice.to_dict(),
                "basis": configuration.crystal.basis.to_dict(),
            }
        )
        elements = configuration.crystal.basis.elements.values
        new_material.basis.remove_atoms_by_elements(elements)
        new_material.basis.add_atom(
            element=configuration.element.chemical_element.value,
            coordinate=configuration.coordinate,
        )
        return new_material


class PointDefectBuilder(MergeBuilder):
    _ConfigurationType: Type[PointDefectConfiguration] = PointDefectConfiguration

    def _merge_component_to_material(self, configuration_or_material: Any) -> Material:
        if isinstance(configuration_or_material, PointDefectSiteConfiguration):
            return PointDefectSiteBuilder().get_material(configuration_or_material)
        return super()._merge_component_to_material(configuration_or_material)

    def _update_material_name(self, material: Material, configuration: _ConfigurationType) -> Material:
        host_material = None
        for component in configuration.merge_components:
            if isinstance(component, Material):
                host_material = component
                break

        if host_material:
            defect_type = configuration.type.lower().replace("defectconfiguration", "").replace("configuration", "")
            material.name = f"{host_material.name} with {defect_type} defect"

        return material


class VacancyDefectBuilder(PointDefectBuilder):
    _ConfigurationType = VacancyDefectConfiguration

    def _post_process(self, material: MaterialWithBuildMetadata, post_process_parameters=None) -> Material:
        material = super()._post_process(material, post_process_parameters)
        material.basis.remove_atoms_by_elements("Vac")
        return material


class SubstitutionalDefectBuilder(PointDefectBuilder):
    _ConfigurationType = SubstitutionalDefectConfiguration


class InterstitialDefectBuilder(PointDefectBuilder):
    _ConfigurationType = InterstitialDefectConfiguration
