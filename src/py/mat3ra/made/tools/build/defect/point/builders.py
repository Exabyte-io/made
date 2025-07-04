from typing import Any

from mat3ra.made.material import Material
from mat3ra.made.tools.build import BaseSingleBuilder
from mat3ra.made.tools.build.merge.builders import MergeBuilder
from mat3ra.made.tools.build.merge.configuration import MergeConfiguration
from mat3ra.made.tools.build.defect.point.configuration import (
    PointDefectConfiguration,
    PointDefectSite,
    VacancyDefectConfiguration,
    SubstitutionalDefectConfiguration,
    InterstitialDefectConfiguration,
)


class PointDefectSiteBuilder(BaseSingleBuilder):
    """
    Builder class for creating a material from a PointDefectSite configuration.
    """

    _ConfigurationType = PointDefectSite

    def _generate(self, configuration: PointDefectSite) -> Material:
        """Creates a single-atom material from scratch using the reference lattice."""
        lattice = configuration.crystal.lattice
        basis = {
            "elements": [{"id": 1, "value": configuration.element.chemical_element}],
            "coordinates": [{"id": 1, "value": configuration.coordinate}],
        }
        return Material.create({"lattice": lattice.to_dict(), "basis": basis, "name": "Defect Site"})


class PointDefectBuilder(MergeBuilder):
    """
    Builder class for creating point defects by merging materials.
    Based on MergeBuilder, similar to how SlabBuilder is based on Stack2ComponentsBuilder.
    """

    _ConfigurationType = PointDefectConfiguration

    def _configuration_to_material(self, configuration_or_material: Any) -> Material:
        if isinstance(configuration_or_material, PointDefectSite):
            return PointDefectSiteBuilder().get_material(configuration_or_material)
        return super()._configuration_to_material(configuration_or_material)

    def _post_process(self, material: Material, configuration: _ConfigurationType) -> Material:
        if isinstance(configuration, VacancyDefectConfiguration):
            material.basis.remove_atoms_by_values("Vac")
        return material

    def get_material(self, configuration: _ConfigurationType) -> Material:
        self.last_configuration = configuration
        merged_material = super().get_material(configuration)
        processed_material = self._post_process(merged_material, configuration)
        return self._update_material_name(processed_material, configuration)

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


class SubstitutionalDefectBuilder(PointDefectBuilder):
    _ConfigurationType = SubstitutionalDefectConfiguration


class InterstitialDefectBuilder(PointDefectBuilder):
    _ConfigurationType = InterstitialDefectConfiguration
