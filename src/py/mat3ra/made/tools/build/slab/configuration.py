from typing import List, Union

import numpy as np
from mat3ra.esse.models.materials_category_components.entities.auxiliary.two_dimensional.miller_indices import (
    MillerIndicesSchema,
)
from mat3ra.esse.models.materials_category_components.entities.reusable.two_dimensional.crystal_lattice_planes import (
    CrystalLatticePlanesSchema,
)
from pydantic import BaseModel

from mat3ra.made.material import Material
from .termination import Termination
from .utils import generate_miller_supercell_matrix, calculate_rotation_matrix, get_terminations, choose_termination
from ..stack.configuration import StackConfiguration
from ..vacuum.configuration import VacuumConfiguration
from ...operations.core.unary import supercell, edit_cell, orient_cell


class MillerSupercell(BaseModel):
    miller_indices: MillerIndicesSchema


class ConventionalCellConfiguration(BaseModel):
    crystal: Material


class CrystalLatticePlanesConfiguration(MillerSupercell, CrystalLatticePlanesSchema):
    crystal: Material

    @property
    def miller_supercell(self) -> List[List[int]]:
        return generate_miller_supercell_matrix(crystal=self.crystal, miller_indices=self.miller_indices)

    @property
    def rotational_matrix(self) -> List[List[float]]:
        miller_supercell_material = supercell(self.crystal, self.miller_supercell)
        return calculate_rotation_matrix(self.crystal, miller_supercell_material)  # to reorient

    @property
    def terminations(self):
        return get_terminations(self.crystal, self.miller_indices)


class AtomicLayersUnique(CrystalLatticePlanesConfiguration):
    @property
    def surface_supercell(self):
        return supercell(self.crystal, self.miller_supercell)

    @property
    def surface_supercell_rotated(self):
        # Rotate the surface supercell to have the Miller indices in the XY plane
        return orient_cell(self.surface_supercell, self.rotational_matrix)

    @property
    def orthogonal_surface_supercell(self):
        supercell_mat = self.surface_supercell_rotated
        current_vectors = supercell_mat.lattice.vector_arrays

        new_vectors = current_vectors.copy()
        new_material = edit_cell(supercell_mat, new_vectors)

        return new_material

    def get_translation_vector(self, termination: Termination) -> List[float]:
        # Implement logic to calculate translation vector based on termination
        return [0.0, 0.0, 0.0]


class AtomicLayersUniqueRepeatedConfiguration(AtomicLayersUnique):
    termination_top: Termination
    number_of_repetitions: int = 1


class SlabConfiguration(StackConfiguration):
    stack_components: List[
        Union[Material, AtomicLayersUnique, AtomicLayersUniqueRepeatedConfiguration, VacuumConfiguration]
    ]
    xy_supercell_matrix: List[List[int]] = [[1, 0], [0, 1]]
