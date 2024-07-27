from typing import List

from mat3ra.made.material import Material
from mat3ra.made.tools.build import BaseBuilder

from .configuration import DeformationConfiguration
from ...utils import solve_sine_wave_x_prime


class DeformationBuilder(BaseBuilder):
    _ConfigurationType: type(DeformationConfiguration) = DeformationConfiguration  # type: ignore
    _GeneratedItemType: Material = Material

    @staticmethod
    def deform_slab(configuration):
        new_material = configuration.slab.clone()
        # new_material.to_cartesian()
        new_coordinates = []
        for coord in new_material.basis.coordinates.values:
            perturbed_coord = configuration.deformation_function[0](coord)
            new_coordinates.append(perturbed_coord)
        new_basis = new_material.basis.copy()
        new_basis.coordinates.values = new_coordinates
        new_basis.to_crystal()
        new_material.basis = new_basis
        return new_material

    def _generate(self, configuration: _ConfigurationType) -> List[_GeneratedItemType]:
        """Generate materials with applied deformation based on the given configuration."""
        new_material = self.deform_slab(configuration)
        return [new_material]

    def _update_material_name(self, material: Material, configuration: _ConfigurationType) -> Material:
        deformation_details = f"Deformation: {configuration.deformation_function[0].__name__}"
        material.name = f"{material.name} ({deformation_details})"
        return material


class ContinuousDeformationBuilder(DeformationBuilder):

    def deform_slab_continuously(self, configuration):
        new_material = configuration.slab.clone()
        new_material.to_cartesian()
        new_coordinates = []
        for coord in new_material.basis.coordinates.values:
            x_prime = solve_sine_wave_x_prime(
                coord[0],
                configuration.deformation_function[1]["amplitude"],
                configuration.deformation_function[1]["wavelength"],
                configuration.deformation_function[1]["phase"],
            )
            perturbed_coord = configuration.deformation_function[0]([x_prime, coord[1], coord[2]])
            new_coordinates.append(perturbed_coord)

        new_basis = new_material.basis.copy()
        new_basis.coordinates.values = new_coordinates
        new_basis.to_crystal()
        new_material.basis = new_basis

        return new_material

    def _generate(
        self, configuration: DeformationBuilder._ConfigurationType
    ) -> List[DeformationBuilder._GeneratedItemType]:
        """Generate materials with applied continuous deformation based on the given configuration."""
        new_material = self.deform_slab_continuously(configuration)
        return [new_material]
