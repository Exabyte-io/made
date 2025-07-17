from typing import List, Optional, Any

from mat3ra.made.tools.build import BaseBuilder, MaterialWithBuildMetadata, BaseSingleBuilder
from .configuration import PerturbationConfiguration
from .parameters import PerturbationBuildParameters
from ...modify import wrap_to_unit_cell, translate_to_z_level
from ...operations.core.unary import perturb


class PerturbationBuilder(BaseSingleBuilder):
    _ConfigurationType: type(PerturbationConfiguration) = PerturbationConfiguration  # type: ignore
    _BuildParametersType: type(PerturbationBuildParameters) = PerturbationBuildParameters
    _DefaultBuildParameters: PerturbationBuildParameters = PerturbationBuildParameters()
    _PostProcessParametersType: Any = None

    def _generate(self, configuration: PerturbationConfiguration) -> MaterialWithBuildMetadata:
        new_material = translate_to_z_level(configuration.material, "center")
        new_material = perturb(new_material, configuration.perturbation_function_holder)
        new_material = wrap_to_unit_cell(new_material)
        return new_material

    def _update_material_name(
        self, material: MaterialWithBuildMetadata, configuration: _ConfigurationType
    ) -> MaterialWithBuildMetadata:
        perturbation_details = f"Perturbation: {configuration.perturbation_function_holder.function_str}"
        material.name = f"{material.name} ({perturbation_details})"
        return material


#
# class SlabPerturbationBuilder(PerturbationBuilder):
#     def create_perturbed_slab(self, configuration: PerturbationConfiguration) -> MaterialWithBuildMetadata:
#         new_material = self._prepare_material(configuration)
#         new_coordinates = [
#             configuration.perturbation_function_holder.apply_perturbation(coord)
#             for coord in new_material.basis.coordinates.values
#         ]
#         new_material.set_coordinates(new_coordinates)
#         return new_material
#
#
# class DistancePreservingSlabPerturbationBuilder(SlabPerturbationBuilder):
#     def create_perturbed_slab(self, configuration: PerturbationConfiguration) -> MaterialWithBuildMetadata:
#         new_material = self._prepare_material(configuration)
#         new_coordinates = [
#             configuration.perturbation_function_holder.transform_coordinates(coord)
#             for coord in new_material.basis.coordinates.values
#         ]
#         new_coordinates = [
#             configuration.perturbation_function_holder.apply_perturbation(coord) for coord in new_coordinates
#         ]
#         new_material.set_coordinates(new_coordinates)
#         return new_material
#
#
# class CellMatchingDistancePreservingSlabPerturbationBuilder(DistancePreservingSlabPerturbationBuilder):
#     def _transform_lattice_vectors(self, configuration: PerturbationConfiguration) -> List[List[float]]:
#         cell_vectors = configuration.material.lattice.vector_arrays
#         return [configuration.perturbation_function_holder.transform_coordinates(coord) for coord in cell_vectors]
#
#     def create_perturbed_slab(self, configuration: PerturbationConfiguration) -> MaterialWithBuildMetadata:
#         new_material = super().create_perturbed_slab(configuration)
#         new_lattice_vectors = self._transform_lattice_vectors(configuration)
#         new_lattice = new_material.lattice.model_copy()
#         new_lattice = new_lattice.from_vectors_array(new_lattice_vectors)
#         new_material.lattice = new_lattice
#         return new_material
