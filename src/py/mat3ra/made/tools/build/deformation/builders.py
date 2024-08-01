from typing import List

from mat3ra.made.material import Material
from mat3ra.made.tools.build import BaseBuilder

from .configuration import DeformationConfiguration
from ...modify import wrap_material


class DeformationBuilder(BaseBuilder):
    _ConfigurationType: type(DeformationConfiguration) = DeformationConfiguration  # type: ignore
    _GeneratedItemType: Material = Material


class SlabDeformationBuilder(DeformationBuilder):
    @staticmethod
    def deform_slab(configuration):
        new_material = configuration.slab.clone()
        new_material.to_cartesian()
        new_coordinates = []
        for coord in new_material.basis.coordinates.values:
            perturbed_coord = configuration.deformation_function[0](coord)
            new_coordinates.append(perturbed_coord)
        new_basis = new_material.basis.copy()
        new_basis.coordinates.values = new_coordinates
        new_basis.to_crystal()
        new_material.basis = new_basis
        return new_material

    def _generate(
        self, configuration: DeformationBuilder._ConfigurationType
    ) -> List[DeformationBuilder._GeneratedItemType]:
        """Generate materials with applied deformation based on the given configuration."""
        new_material = self.deform_slab(configuration)
        return [new_material]

    def _update_material_name(
        self, material: Material, configuration: DeformationBuilder._ConfigurationType
    ) -> Material:
        deformation_details = f"Deformation: {configuration.deformation_function[0].__name__}"
        material.name = f"{material.name} ({deformation_details})"
        return material


class DistancePreservingSlabDeformationBuilder(DeformationBuilder):
    def deform_slab_isometrically(self, configuration):
        new_material = configuration.slab.clone()
        if configuration.use_cartesian_coordinates:
            new_material.to_cartesian()
        new_coordinates = []

        deformation_function, deformation_json = configuration.deformation_function
        for coord in new_material.basis.coordinates.values:
            perturbed_coord = deformation_function(coord)
            new_coordinates.append(perturbed_coord)

        new_basis = new_material.basis.copy()
        new_basis.coordinates.values = new_coordinates
        new_basis.to_crystal()
        new_material.basis = new_basis

        new_material = wrap_material(new_material)
        return new_material

    def _generate(
        self, configuration: DeformationBuilder._ConfigurationType
    ) -> List[DeformationBuilder._GeneratedItemType]:
        """Generate materials with applied continuous deformation based on the given configuration."""
        new_material = self.deform_slab_isometrically(configuration)
        return [new_material]
