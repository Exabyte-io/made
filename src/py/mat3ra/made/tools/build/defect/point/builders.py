from typing import Any, Type

from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from mat3ra.esse.models.materials_category_components.entities.core.zero_dimensional.atom import (
    AtomSchema,
)

from mat3ra.made.material import Material
from mat3ra.made.tools.build import BaseSingleBuilder, MaterialWithBuildMetadata
from mat3ra.made.tools.build.defect.point.configuration import (
    PointDefectConfiguration,
    PointDefectSiteConfiguration,
    VacancyDefectConfiguration,
    SubstitutionalDefectConfiguration,
    InterstitialDefectConfiguration,
)
from mat3ra.made.tools.build.merge.builders import MergeBuilder
from mat3ra.made.tools.build.merge.configuration import MergeConfiguration
from mat3ra.made.tools.build.vacuum.builders import VacuumBuilder
from mat3ra.made.tools.build.vacuum.configuration import VacuumConfiguration


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
    """
    Builder class for creating point defects by merging materials.
    Based on MergeBuilder, similar to how SlabBuilder is based on StackNComponentsBuilder.
    """

    _ConfigurationType: Type[PointDefectConfiguration] = PointDefectConfiguration

    def _configuration_to_material(self, configuration_or_material: Any) -> Material:
        if isinstance(configuration_or_material, PointDefectSiteConfiguration):
            return PointDefectSiteBuilder().get_material(configuration_or_material)
        return super()._configuration_to_material(configuration_or_material)

    def _post_process(self, material: Material, configuration: MergeBuilder._ConfigurationType) -> Material:
        if isinstance(configuration, VacancyDefectConfiguration):
            material.basis.remove_atoms_by_elements("Vac")
        return material

    def get_material(self, configuration: MergeBuilder._ConfigurationType) -> Material:
        merged_material = super().get_material(configuration)
        processed_material = self._post_process(merged_material, configuration)
        return self._update_material_name(processed_material, configuration)

    def _update_material_name(self, material: Material, configuration: MergeBuilder._ConfigurationType) -> Material:
        host_material = None
        for component in configuration.merge_components:
            if isinstance(component, Material):
                host_material = component
                break

        if host_material:
            defect_type = configuration.type.lower().replace("defectconfiguration", "").replace("configuration", "")
            material.name = f"{host_material.name} with {defect_type} defect"

        return material

    def _generate(self, config: MergeBuilder._ConfigurationType) -> MaterialWithBuildMetadata:
        materials = []
        site_builder = PointDefectSiteBuilder()
        for component in config.merge_components:
            if isinstance(component, PointDefectSiteConfiguration):
                materials.append(site_builder.get_material(component))
            elif isinstance(component, Material):
                materials.append(component)

        merge_config = MergeConfiguration(
            merge_components=materials,
            merge_method=config.merge_method.value,
        )

        return super()._generate(merge_config)


class VacancyDefectBuilder(PointDefectBuilder):
    _ConfigurationType = VacancyDefectConfiguration


class SubstitutionalDefectBuilder(PointDefectBuilder):
    _ConfigurationType = SubstitutionalDefectConfiguration


class InterstitialDefectBuilder(PointDefectBuilder):
    _ConfigurationType = InterstitialDefectConfiguration
