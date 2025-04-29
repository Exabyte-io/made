from typing import List, Optional, Any

from mat3ra.made.material import Material
from mat3ra.made.tools.build import BaseBuilder

from .configuration import PerturbationConfiguration
from ...modify import wrap_to_unit_cell, translate_to_z_level


class PerturbationBuilder(BaseBuilder):
    _ConfigurationType: type(PerturbationConfiguration) = PerturbationConfiguration  # type: ignore
    _GeneratedItemType: Material = Material
    _PostProcessParametersType: Any = None

    @staticmethod
    def _prepare_material(configuration: _ConfigurationType) -> _GeneratedItemType:
        new_material = configuration.material.clone()
        new_material = translate_to_z_level(new_material, "center")
        if configuration.use_cartesian_coordinates:
            new_material.to_cartesian()
        return new_material

    def _generate(self, configuration: _ConfigurationType) -> List[_GeneratedItemType]:
        """Generate materials with applied continuous perturbation based on the given configuration."""
        new_material = self.create_perturbed_slab(configuration)
        return [new_material]

    def _post_process(
        self,
        items: List[_GeneratedItemType],
        post_process_parameters: Optional[_PostProcessParametersType],
    ) -> List[Material]:
        return [wrap_to_unit_cell(item) for item in items]

    def _update_material_name(self, material: Material, configuration: _ConfigurationType) -> Material:
        perturbation_details = f"Perturbation: {configuration.perturbation_function_holder.get_json().get('type')}"
        material.name = f"{material.name} ({perturbation_details})"
        return material


class SlabPerturbationBuilder(PerturbationBuilder):
    def create_perturbed_slab(self, configuration: PerturbationConfiguration) -> Material:
        new_material = self._prepare_material(configuration)
        new_coordinates = [
            configuration.perturbation_function_holder.apply_perturbation(coord)
            for coord in new_material.basis.coordinates.values
        ]
        new_material.set_coordinates(new_coordinates)
        return new_material


class DistancePreservingSlabPerturbationBuilder(SlabPerturbationBuilder):
    def create_perturbed_slab(self, configuration: PerturbationConfiguration) -> Material:
        new_material = self._prepare_material(configuration)
        new_coordinates = [
            configuration.perturbation_function_holder.transform_coordinates(coord)
            for coord in new_material.basis.coordinates.values
        ]
        new_coordinates = [
            configuration.perturbation_function_holder.apply_perturbation(coord) for coord in new_coordinates
        ]
        new_material.set_coordinates(new_coordinates)
        return new_material


class CellMatchingDistancePreservingSlabPerturbationBuilder(DistancePreservingSlabPerturbationBuilder):
    def _transform_lattice_vectors(self, configuration: PerturbationConfiguration) -> List[List[float]]:
        cell_vectors = configuration.material.lattice.vector_arrays
        return [configuration.perturbation_function_holder.transform_coordinates(coord) for coord in cell_vectors]

    def create_perturbed_slab(self, configuration: PerturbationConfiguration) -> Material:
        new_material = super().create_perturbed_slab(configuration)
        new_lattice_vectors = self._transform_lattice_vectors(configuration)
        new_lattice = new_material.lattice.model_copy()
        new_lattice = new_lattice.from_vectors_array(new_lattice_vectors)
        new_material.lattice = new_lattice
        return new_material
