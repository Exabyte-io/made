from typing import Optional

from mat3ra.esse.models.materials_category.defective_structures.two_dimensional.base_configuration import (
    AtomicLayersUniqueRepeatedSchema,
)

from mat3ra.made.tools.build.slab.crystal_lattice_planes.configuration import CrystalLatticePlanesConfiguration
from mat3ra.made.tools.build.slab.entities import Termination


class AtomicLayersUniqueRepeatedConfiguration(CrystalLatticePlanesConfiguration, AtomicLayersUniqueRepeatedSchema):
    termination_top: Optional[Termination]
    number_of_repetitions: int
