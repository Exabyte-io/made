from typing import Optional

from mat3ra.esse.models.materials_category_components.entities.reusable.two_dimensional.crystal_lattice_planes import (
    CrystalLatticePlanesSchema,
)

from mat3ra.made.material import Material
from mat3ra.made.tools.build.slab.entities import Termination
from ... import BaseConfigurationPydantic


class CrystalLatticePlanesConfiguration(CrystalLatticePlanesSchema, BaseConfigurationPydantic):
    crystal: Material

    @property
    def in_plane_vectors(self):
        # Two vectors in the plane of the Miller indices
        return self.crystal.lattice.vector_arrays[:2, :2]


class AtomicLayersUnique(CrystalLatticePlanesConfiguration):
    pass


class AtomicLayersUniqueRepeatedConfiguration(AtomicLayersUnique):
    termination_top: Termination
    termination_bottom: Optional[Termination] = None
    number_of_repetitions: int = 1 