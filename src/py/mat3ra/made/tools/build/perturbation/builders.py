from typing import List

from mat3ra.made.material import Material
from mat3ra.made.tools.build import BaseBuilder

from .configuration import PerturbationConfiguration
from ...modify import wrap_material
from ...utils import PerturbationFunctionHolder


class PerturbationBuilder(BaseBuilder):
    _ConfigurationType: type(PerturbationConfiguration) = PerturbationConfiguration  # type: ignore
    _GeneratedItemType: Material = Material

    @staticmethod
    def _prepare_material(configuration):
        new_material = configuration.material.clone()
        if configuration.use_cartesian_coordinates:
            new_material.to_cartesian()
        return new_material

    @staticmethod
    def _set_new_coordinates(new_material, new_coordinates):
        new_basis = new_material.basis.copy()
        new_basis.coordinates.values = new_coordinates
        new_basis.to_crystal()
        new_material.basis = new_basis
        return new_material

    def _generate(self, configuration: _ConfigurationType) -> List[_GeneratedItemType]:
        """Generate materials with applied continuous perturbation based on the given configuration."""
        new_material = self.create_perturbed_slab(configuration)
        return [new_material]

    def _update_material_name(self, material: Material, configuration: _ConfigurationType) -> Material:
        perturbation_details = f"Perturbation: {configuration.perturbation_function[1].get('type')}"
        material.name = f"{material.name} ({perturbation_details})"
        return material


class SlabPerturbationBuilder(PerturbationBuilder):
    def create_perturbed_slab(self, configuration):
        new_material = self._prepare_material(configuration)
        new_coordinates = []
        perturbation_function, _ = configuration.perturbation_function
        for coord in new_material.basis.coordinates.values:
            perturbed_coord = perturbation_function(coord)
            new_coordinates.append(perturbed_coord)
        new_material = self._set_new_coordinates(new_material, new_coordinates)
        return new_material


class DistancePreservingSlabPerturbationBuilder(PerturbationBuilder):
    def create_perturbed_slab(self, configuration):
        new_material = self._prepare_material(configuration)
        new_coordinates = []
        perturbation_function, perturbation_json = configuration.perturbation_function
        coord_transformation_function = PerturbationFunctionHolder.get_coord_transformation(perturbation_json)
        for coord in new_material.basis.coordinates.values:
            transformed_coord = coord_transformation_function(coord)
            perturbed_coord = perturbation_function(transformed_coord)
            new_coordinates.append(perturbed_coord)
        new_material = self._set_new_coordinates(new_material, new_coordinates)
        new_material = wrap_material(new_material)
        return new_material
