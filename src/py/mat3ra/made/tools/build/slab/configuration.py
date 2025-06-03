from typing import List, Union

from mat3ra.esse.models.material.primitive.two_dimensional.miller_indices import MillerIndicesSchema
from mat3ra.esse.models.materials_category.pristine_structures.two_dimensional.slab import (
    VacuumConfigurationSchema,
)
from pydantic import BaseModel

from mat3ra.made.material import Material
from .termination import Termination
from .. import BaseBuilder
from ...convert import to_pymatgen
from ...operations.core.unary import supercell, stack
from ...third_party import PymatgenSlab
from ...third_party import PymatgenSlabGenerator, label_pymatgen_slab_termination


def generate_pymatgen_slabs(
    crystal: Material,
    miller_indices: Union[MillerIndicesSchema, List[int]] = (0, 0, 1),
    min_slab_size=1,
    min_vacuum_size=0,
    in_unit_planes: bool = True,
    make_primitive: bool = False,
    symmetrize: bool = False,
) -> List[PymatgenSlab]:  # type: ignore
    generator = PymatgenSlabGenerator(
        initial_structure=to_pymatgen(crystal),
        miller_index=miller_indices,
        min_slab_size=min_slab_size,
        min_vacuum_size=min_vacuum_size,
        in_unit_planes=in_unit_planes,
        primitive=make_primitive,
    )
    raw_slabs = generator.get_slabs(
        # We need to preserve symmetric slabs for different terminations at the surface
        symmetrize=symmetrize
    )

    return raw_slabs


class MillerSupercell(BaseModel):
    miller_indices: MillerIndicesSchema

    @property
    def miller_supercell(self):
        return ...

        # use pymatgen


def calculate_rotation_matrix(crystal, miller_supercell_material):

    # Implement logic to calculate the rotation matrix based on crystal and miller_supercell_material
    return [[1, 0, 0], [0, 1, 0], [0, 0, 1]]  # Identity matrix as a placeholder


def get_terminations(crystal: Material, miller_indices: Union[MillerIndicesSchema, List[int]]) -> List[Termination]:
    return [
        Termination.from_string(label_pymatgen_slab_termination(slab))
        for slab in generate_pymatgen_slabs(crystal, miller_indices)
    ]


def choose_termination(terminations: List[Termination], stoichiometry: str) -> Termination:
    # choose a termination by stochiometry or symmetry provided
    for termination in terminations:
        if str(termination.chemical_elements) == name:
            return termination
    return terminations[0] if terminations else None


def create_vacuum_material(reference: Material, vacuum: "VacuumConfiguration") -> Material:
    """
    Create a vacuum material based on a reference material and vacuum configuration.

    Args:
        reference: Reference material to base the vacuum lattice on.
        vacuum: Vacuum configuration with size and direction.

    Returns:
        Material: Vacuum material with empty basis.
    """

    a_vector, b_vector = reference.lattice.vector_arrays[:2]
    vacuum_lattice = reference.lattice.from_vectors_array(
        [a_vector, b_vector, [0, 0, vacuum.size]], reference.lattice.units, reference.lattice.type
    )
    return Material.create(
        {
            "name": "Vacuum",
            "lattice": vacuum_lattice.to_dict(),
            "basis": {"elements": [], "coordinates": []},
        }
    )


class OrientedCrystal(MillerSupercell):
    crystal: Material

    @property
    def rotational_matrix(self) -> List[List[float]]:
        miller_supercell_material = supercell(self.crystal, self.miller_supercell)
        return calculate_rotation_matrix(self.crystal, miller_supercell_material)  # to reorient

    @property
    def terminations(self):
        return get_terminations(self.crystal, self.miller_indices)


def calculate_orthogonal_lattice(crystal, miller_supercell):
    # Implement logic to calculate the orthogonal lattice based on crystal and miller_supercell
    return supercell(crystal, miller_supercell)  # Placeholder implementation


class AtomicLayersUnique(OrientedCrystal):
    @property
    def orthogonal_c_cell(self):
        return calculate_orthogonal_lattice(self.crystal, self.miller_supercell)

    def get_translation_vector(self, termination: Termination) -> List[float]:
        # Implement logic to calculate translation vector based on termination
        return [0.0, 0.0, 0.0]


class AtomicLayersUniqueRepeated(AtomicLayersUnique):
    number_of_repetitions: int

    @property
    def repeated_layers(self) -> Material:
        """
        Returns the repeated cell based on the number of repetitions.
        """
        return supercell(self.orthogonal_c_cell, [[1, 0, 0], [0, 1, 0], [0, 0, self.number_of_repetitions]])


class VacuumConfiguration(VacuumConfigurationSchema):
    size: float

    @property
    def vacuum_layer(self, reference: Material) -> Material:
        return create_vacuum_material(reference, self.size)


class SlabConfiguration(BaseModel):
    atomic_layers: AtomicLayersUniqueRepeated
    vacuum_configuration: VacuumConfigurationSchema

    def from_parameters(
        crystal: Material,
        miller_indices: Union[MillerIndicesSchema, List[int]] = (0, 0, 1),
        use_conventional_cell: bool = True,
        termination: Termination = None,
        number_of_layers: int = 1,
        vacuum: float = 10.0,
        xy_supercell_matrix: List[List[int]] = None,
    ) -> "SlabConfiguration":
        atomic_layers = AtomicLayersUniqueRepeated(
            crystal=crystal,
            miller_indices=miller_indices,
            use_conventional_cell=use_conventional_cell,
            number_of_repetitions=number_of_layers,
        )
        vacuum_configuration = VacuumConfiguration(size=vacuum)
        return SlabConfiguration(atomic_layers=atomic_layers, vacuum_configuration=vacuum_configuration)


class SlabBuilder(BaseBuilder):

    def generate(self, configuration: SlabConfiguration) -> Material:
        repeated_layers = configuration.atomic_layers.repeated_layers
        vacuum_layer = configuration.vacuum_configuration.vacuum_layer
        stacked_materials = stack([repeated_layers, vacuum_layer], "z")
        supercell_slab = supercell(stacked_materials, configuration.xy_supercell_matrix)
        return supercell_slab

    def get_material(self, configuration: SlabConfiguration) -> Material:
        return self.generate(configuration)
