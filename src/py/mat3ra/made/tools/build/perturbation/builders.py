from typing import List

from mat3ra.made.material import Material
from mat3ra.made.tools.build import BaseBuilder

from .configuration import PerturbationConfiguration
from ...modify import wrap_material


class PerturbationBuilder(BaseBuilder):
    _ConfigurationType: type(PerturbationConfiguration) = PerturbationConfiguration  # type: ignore
    _GeneratedItemType: Material = Material


class SlabPerturbationBuilder(PerturbationBuilder):
    @staticmethod
    def deform_slab(configuration):
        new_material = configuration.material.clone()
        new_material.to_cartesian()
        new_coordinates = []
        for coord in new_material.basis.coordinates.values:
            perturbed_coord = configuration.perturbation_function[0](coord)
            new_coordinates.append(perturbed_coord)
        new_basis = new_material.basis.copy()
        new_basis.coordinates.values = new_coordinates
        new_basis.to_crystal()
        new_material.basis = new_basis
        return new_material

    def _generate(
        self, configuration: PerturbationBuilder._ConfigurationType
    ) -> List[PerturbationBuilder._GeneratedItemType]:
        """Generate materials with applied perturbation based on the given configuration."""
        new_material = self.deform_slab(configuration)
        return [new_material]

    def _update_material_name(
        self, material: Material, configuration: PerturbationBuilder._ConfigurationType
    ) -> Material:
        perturbation_details = f"Perturbation: {configuration.perturbation_function[0].__name__}"
        material.name = f"{material.name} ({perturbation_details})"
        return material


class DistancePreservingSlabPerturbationBuilder(PerturbationBuilder):
    def deform_slab_isometrically(self, configuration):
        new_material = configuration.material.clone()
        if configuration.use_cartesian_coordinates:
            new_material.to_cartesian()
        new_coordinates = []

        perturbation_function, perturbation_json = configuration.perturbation_function
        for coord in new_material.basis.coordinates.values:
            perturbed_coord = perturbation_function(coord)
            new_coordinates.append(perturbed_coord)

        new_basis = new_material.basis.copy()
        new_basis.coordinates.values = new_coordinates
        new_basis.to_crystal()
        new_material.basis = new_basis

        new_material = wrap_material(new_material)
        return new_material

    def _generate(
        self, configuration: PerturbationBuilder._ConfigurationType
    ) -> List[PerturbationBuilder._GeneratedItemType]:
        """Generate materials with applied continuous perturbation based on the given configuration."""
        new_material = self.deform_slab_isometrically(configuration)
        return [new_material]
