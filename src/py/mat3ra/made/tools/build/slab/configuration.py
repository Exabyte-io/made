from typing import List, Union

from mat3ra.esse.models.materials_category_components.entities.auxiliary.two_dimensional.miller_indices import (
    MillerIndicesSchema,
)
from mat3ra.esse.models.materials_category_components.entities.reusable.two_dimensional.crystal_lattice_planes import (
    CrystalLatticePlanesSchema,
)
from pydantic import BaseModel

from mat3ra.made.material import Material
from .termination import Termination
from .utils import (
    generate_miller_supercell_matrix,
    get_terminations,
)
from ..stack.configuration import StackConfiguration
from ..vacuum.configuration import VacuumConfiguration
from ...operations.core.unary import supercell


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
    def terminations(self):
        return get_terminations(self.crystal, self.miller_indices)


class AtomicLayersUnique(CrystalLatticePlanesConfiguration):
    @property
    def surface_supercell(self):
        return supercell(self.crystal, self.miller_supercell)

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
    use_orthogonal_c: bool = True
