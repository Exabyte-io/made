from typing import Any

from mat3ra.made.tools.build import MaterialWithBuildMetadata, BaseSingleBuilder
from .configuration import PerturbationConfiguration
from .parameters import PerturbationBuildParameters
from ...modify import wrap_to_unit_cell, translate_to_z_level
from ...operations.core.unary import perturb, edit_cell


class PerturbationBuilder(BaseSingleBuilder):
    _ConfigurationType: type(PerturbationConfiguration) = PerturbationConfiguration  # type: ignore
    _BuildParametersType: type(PerturbationBuildParameters) = PerturbationBuildParameters
    _DefaultBuildParameters: PerturbationBuildParameters = PerturbationBuildParameters()
    _PostProcessParametersType: Any = None

    def _generate(self, configuration: PerturbationConfiguration) -> MaterialWithBuildMetadata:
        new_material = translate_to_z_level(configuration.material, "center").clone()
        new_material = perturb(
            new_material,
            configuration.perturbation_function_holder,
            configuration.use_cartesian_coordinates,
        )
        new_material = wrap_to_unit_cell(new_material)
        return new_material

    def _update_material_name(
        self, material: MaterialWithBuildMetadata, configuration: _ConfigurationType
    ) -> MaterialWithBuildMetadata:
        perturbation_details = f"Perturbation: {configuration.perturbation_function_holder.function_str}"
        material.name = f"{material.name} ({perturbation_details})"
        return material


class IsometricPerturbationBuilder(PerturbationBuilder):
    def _generate(self, configuration: PerturbationConfiguration) -> MaterialWithBuildMetadata:
        new_material = configuration.material.clone()
        renormalized_coordinates = [
            configuration.perturbation_function_holder.normalize_coordinates(coord)
            for coord in new_material.basis.coordinates.values
        ]
        new_material.set_coordinates(renormalized_coordinates)

        configuration.material = new_material

        new_material = super()._generate(configuration)

        new_lattice_vectors = [
            configuration.perturbation_function_holder.normalize_coordinates(vector)
            for vector in new_material.lattice.vector_arrays
        ]
        renormalized_material = edit_cell(new_material, new_lattice_vectors)
        return renormalized_material
