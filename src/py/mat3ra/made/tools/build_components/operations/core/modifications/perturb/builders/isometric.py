from .......operations.core.unary import edit_cell
from ...... import MaterialWithBuildMetadata
from ..configuration import PerturbationConfiguration
from .base import PerturbationBuilder


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
