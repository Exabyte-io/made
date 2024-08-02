from typing import List, Optional

from mat3ra.made.material import Material
from mat3ra.made.tools.build import BaseBuilder

from .configuration import PerturbationConfiguration
from ...modify import wrap_material, translate_to_z_level
from ...utils import PerturbationFunctionHolder


class PerturbationBuilder(BaseBuilder):
    _ConfigurationType: type(PerturbationConfiguration) = PerturbationConfiguration  # type: ignore
    _GeneratedItemType: Material = Material
    _PostProcessParametersType = None

    @staticmethod
    def _prepare_material(configuration):
        new_material = configuration.material.clone()
        new_material = translate_to_z_level(new_material, "center")
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

    def _post_process(
        self, items: List[_GeneratedItemType], post_process_parameters: Optional[_PostProcessParametersType]
    ) -> List[Material]:
        return [wrap_material(item) for item in items]

    def _update_material_name(self, material: Material, configuration: _ConfigurationType) -> Material:
        perturbation_details = f"Perturbation: {configuration.perturbation_function[1].get('type')}"
        material.name = f"{material.name} ({perturbation_details})"
        return material


class SlabPerturbationBuilder(PerturbationBuilder):
    def create_perturbed_slab(self, configuration):
        new_material = self._prepare_material(configuration)
        perturbation_function, _ = configuration.perturbation_function
        new_coordinates = [perturbation_function(coord) for coord in new_material.basis.coordinates.values]
        new_material = self._set_new_coordinates(new_material, new_coordinates)
        return new_material


class DistancePreservingSlabPerturbationBuilder(SlabPerturbationBuilder):
    def create_perturbed_slab(self, configuration):
        new_material = self._prepare_material(configuration)
        perturbation_function, perturbation_json = configuration.perturbation_function
        coord_transformation_function = PerturbationFunctionHolder.get_coord_transformation(perturbation_json)
        new_coordinates = [
            perturbation_function(coord_transformation_function(coord))
            for coord in new_material.basis.coordinates.values
        ]
        new_material = self._set_new_coordinates(new_material, new_coordinates)
        return new_material


class CellMatchingDistancePreservingSlabPerturbationBuilder(DistancePreservingSlabPerturbationBuilder):
    def _transform_cell_vectors(self, configuration: PerturbationConfiguration) -> List[List[float]]:
        perturbation_function, perturbation_json = configuration.perturbation_function
        coord_transformation_function = PerturbationFunctionHolder.get_coord_transformation(perturbation_json)
        cell_vectors = configuration.material.basis.cell.vectors_as_nested_array
        return [perturbation_function(coord_transformation_function(coord)) for coord in cell_vectors]

    def create_perturbed_slab(self, configuration: PerturbationConfiguration):
        new_material = super().create_perturbed_slab(configuration)
        new_lattice_vectors = self._transform_cell_vectors(configuration)
        new_lattice = new_material.lattice.copy()
        new_lattice = new_lattice.from_nested_array(new_lattice_vectors)
        new_material.lattice = new_lattice

        new_basis = new_material.basis.copy()
        new_basis.to_cartesian()
        new_basis.cell = new_basis.cell.from_nested_array(new_lattice_vectors)
        new_basis.to_crystal()
        new_material.basis = new_basis
        return new_material
