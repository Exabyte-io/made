from typing import List, Union

from mat3ra.code.vector import Vector3D
from mat3ra.esse.models.materials_category_components.entities.auxiliary.two_dimensional.miller_indices import (
    MillerIndicesSchema,
)
from mat3ra.esse.models.materials_category_components.entities.reusable.two_dimensional.crystal_lattice_planes import (
    CrystalLatticePlanesSchema,
)
from pydantic import BaseModel

from mat3ra.made.material import Material
from .termination import Termination
from .. import BaseBuilder
from ..stack.configuration import StackConfiguration
from ..utils import miller_to_supercell_matrix
from ..vacuum.configuration import VacuumConfiguration
from ...convert import to_pymatgen, from_pymatgen
from ...operations.core.unary import supercell, translate
from ...third_party import PymatgenSlab
from ...third_party import PymatgenSlabGenerator, label_pymatgen_slab_termination, PymatgenSpacegroupAnalyzer


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


class ConventionalCellConfiguration(BaseModel):
    crystal: Material


class ConventionalCellBuilder(BaseBuilder):
    def get_material(self, configuration: ConventionalCellConfiguration) -> Material:
        bulk_pmg = to_pymatgen(configuration.crystal)
        conventional_cell = PymatgenSpacegroupAnalyzer(bulk_pmg).get_conventional_standard_structure()
        return Material.create(from_pymatgen(conventional_cell))


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
    @staticmethod
    def calculate_orthogonal_lattice(crystal, miller_supercell):
        return supercell(crystal, miller_supercell)

    def get_material(self, configuration: CrystalLatticePlanesConfiguration) -> Material:
        orthogonal_c_cell = self.calculate_orthogonal_lattice(configuration.crystal, configuration.miller_supercell)
        return supercell(orthogonal_c_cell, configuration.rotational_matrix)


class AtomicLayersUnique(CrystalLatticePlanesConfiguration):
    @property
    def orthogonal_c_cell(self):
        return supercell(self.crystal, self.miller_supercell)

    def get_translation_vector(self, termination: Termination) -> List[float]:
        # Implement logic to calculate translation vector based on termination
        return [0.0, 0.0, 0.0]


class AtomicLayersUniqueRepeatedConfiguration(AtomicLayersUnique):
    termination_top: Termination
    number_of_repetitions: int = 1


class AtomicLayersUniqueRepeatedBuilder(BaseBuilder):
    _ConfigurationType = AtomicLayersUniqueRepeatedConfiguration

    def get_material(self, configuration: _ConfigurationType) -> Material:
        material = configuration.orthogonal_c_cell
        translation_vector: Vector3D = configuration.get_translation_vector(configuration.termination_top)
        material_translated = translate(material, translation_vector)
        material_translated_with_repetitions = supercell(
            material_translated, [[1, 0, 0], [0, 1, 0], [0, 0, configuration.number_of_repetitions]]
        )
        return material_translated_with_repetitions


class SlabConfiguration(StackConfiguration):
    stack_components: List[
        Union[Material, AtomicLayersUnique, AtomicLayersUniqueRepeatedConfiguration, VacuumConfiguration]
    ]
    xy_supercell_matrix: List[List[int]] = [[1, 0], [0, 1]]
