from typing import List, Union, Any

from mat3ra.esse.models.material.primitive.combinations.stack import StackSchema
from mat3ra.esse.models.material.primitive.two_dimensional.miller_indices import MillerIndicesSchema
from mat3ra.esse.models.material.reusable.two_dimensional.crystal_lattice_planes import (
    CrystalLatticePlanesSchema,
)
from mat3ra.esse.models.materials_category.pristine_structures.two_dimensional.slab import (
    SlabConfigurationSchema,
    AxisEnum,
)
from pydantic import BaseModel

from mat3ra.made.material import Material
from .termination import Termination
from .. import BaseBuilder
from ..utils import miller_to_supercell_matrix
from ..vacuum.builders import VacuumBuilder
from ..vacuum.configuration import VacuumConfiguration
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
    # Extract actual values from MillerIndicesSchema if needed
    if isinstance(miller_indices, MillerIndicesSchema):
        miller_values = miller_indices.root
    else:
        miller_values = miller_indices

    generator = PymatgenSlabGenerator(
        initial_structure=to_pymatgen(crystal),
        miller_index=miller_values,
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
        if str(termination.chemical_elements) == stoichiometry:
            return termination
    return terminations[0] if terminations else None


class MillerSupercell(BaseModel):
    miller_indices: MillerIndicesSchema

    @property
    def miller_supercell(self):
        # use pymatgen to generate the supercell based on miller indices
        return miller_to_supercell_matrix(self.miller_indices.root)


class CrystalLatticePlanesConfiguration(MillerSupercell, CrystalLatticePlanesSchema):
    crystal: Material

    @property
    def rotational_matrix(self) -> List[List[float]]:
        miller_supercell_material = supercell(self.crystal, self.miller_supercell)
        return calculate_rotation_matrix(self.crystal, miller_supercell_material)  # to reorient

    @property
    def terminations(self):
        return get_terminations(self.crystal, self.miller_indices)


class CrystalLatticePlanesBuilder(BaseBuilder):
    def calculate_orthogonal_lattice(crystal, miller_supercell):
        # Implement logic to calculate the orthogonal lattice based on crystal and miller_supercell
        return supercell(crystal, miller_supercell)  # Placeholder implementation

    def get_material(self, configuration: CrystalLatticePlanesConfiguration) -> Material:
        orthogonal_c_cell = self.calculate_orthogonal_lattice(configuration.crystal, configuration.miller_supercell)
        return supercell(orthogonal_c_cell, configuration.rotational_matrix)


class AtomicLayersUnique(CrystalLatticePlanesConfiguration):
    @property
    def orthogonal_c_cell(self):
        # return calculate_orthogonal_lattice(self.crystal, self.miller_supercell)
        return supercell(self.crystal, [[1, 0, 0], [0, 1, 0], [0, 0, 1]])  # Placeholder implementation

    def get_translation_vector(self, termination: Termination) -> List[float]:
        # Implement logic to calculate translation vector based on termination
        return [0.0, 0.0, 0.0]


class AtomicLayersUniqueRepeatedConfiguration(AtomicLayersUnique):
    number_of_repetitions: int = 1


class AtomicLayersUniqueRepeatedBuilder(BaseBuilder):
    _ConfigurationType = AtomicLayersUniqueRepeatedConfiguration

    def get_material(self, configuration: _ConfigurationType) -> Material:
        """
        Returns the repeated cell based on the number of repetitions.
        """
        return supercell(
            configuration.orthogonal_c_cell, [[1, 0, 0], [0, 1, 0], [0, 0, configuration.number_of_repetitions]]
        )


class StackConfiguration(StackSchema):

    stack_components: List[Union[AtomicLayersUniqueRepeatedConfiguration, VacuumConfiguration]]

    @property
    def atomic_layers(self) -> AtomicLayersUniqueRepeatedConfiguration:
        return self.stack_components[0]

    @property
    def vacuum_configuration(self) -> VacuumConfiguration:
        return self.stack_components[1]


class SlabConfiguration(SlabConfigurationSchema):
    stack_components: List[Union[AtomicLayersUnique, VacuumConfiguration]]


class StackBuilder2Components(BaseBuilder):

    def configuration_to_material(self, configuration_or_material: Any) -> Material:
        if isinstance(configuration_or_material, Material):
            return configuration_or_material
        if isinstance(configuration_or_material, AtomicLayersUniqueRepeatedConfiguration):
            builder = AtomicLayersUniqueRepeatedBuilder()
        if isinstance(configuration_or_material, VacuumConfiguration):
            builder = VacuumBuilder()
        return builder.get_material(configuration_or_material)

    def generate(self, configuration: StackConfiguration) -> Material:
        first_entity_config = self.configuration.stack_components[0]
        first_material = self.configuration_to_material(first_entity_config)
        second_entity_config = self.configuration.stack_components[1]
        second_material = self.configuration_to_material(second_entity_config)
        # Stack the two materials
        stacked_materials = stack([first_material, second_material], AxisEnum.z)

    def get_material(self, configuration: SlabConfiguration) -> Material:
        return self.generate(configuration)
