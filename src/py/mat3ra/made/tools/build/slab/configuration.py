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
    generate_pymatgen_slabs,
)
from ..stack.configuration import StackConfiguration
from ..vacuum.configuration import VacuumConfiguration


# TODO: should be moved to analyze module
class TerminationAnalyzer:
    """Detects and calculates translation vectors for slab terminations."""

    def __init__(self, crystal: Material, miller_indices: Union[MillerIndicesSchema, List[int]]):
        self.crystal = crystal
        self.miller_indices = miller_indices

    def find_translation_vector(self, target_termination: Termination) -> List[float]:
        """Find translation vector in crystal units to achieve target termination at slab top."""
        from .builders import ConventionalCellBuilder

        conv_config = ConventionalCellConfiguration(crystal=self.crystal)
        conv_crystal = ConventionalCellBuilder().get_material(conv_config)

        slabs = generate_pymatgen_slabs(
            crystal=conv_crystal,
            miller_indices=self.miller_indices,
            min_vacuum_size=5.0,  # Add vacuum for termination detection
        )
        for slab in slabs:
            from ...third_party import label_pymatgen_slab_termination

            slab_termination_str = label_pymatgen_slab_termination(slab)
            slab_termination = Termination.from_string(slab_termination_str)

            if slab_termination == target_termination:
                if hasattr(slab, "shift"):
                    shift = slab.shift
                    return [0.0, 0.0, float(shift)]

        return [0.0, 0.0, 0.0]


class ConventionalCellConfiguration(BaseModel):
    crystal: Material


class CrystalLatticePlanesConfiguration(CrystalLatticePlanesSchema):
    crystal: Material

    @property
    def miller_supercell(self) -> List[List[int]]:
        return generate_miller_supercell_matrix(crystal=self.crystal, miller_indices=self.miller_indices)

    @property
    def terminations(self):
        return get_terminations(self.crystal, self.miller_indices)


class AtomicLayersUnique(CrystalLatticePlanesConfiguration):
    pass


class AtomicLayersUniqueRepeatedConfiguration(AtomicLayersUnique):
    termination_top: Termination
    number_of_repetitions: int = 1


class SlabConfiguration(StackConfiguration):
    stack_components: List[
        Union[AtomicLayersUnique, AtomicLayersUniqueRepeatedConfiguration, VacuumConfiguration]  # No Materials!
    ]

    @property
    def atomic_layers(self):
        return self.stack_components[0]

    @property
    def vacuum_configuration(self) -> VacuumConfiguration:
        return self.stack_components[1]
