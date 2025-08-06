from typing import Optional

from mat3ra.esse.models.materials_category.defective_structures.two_dimensional.base_configuration import (
    AtomicLayersUniqueRepeatedSchema,
)

from ..crystal_lattice_planes.configuration import CrystalLatticePlanesConfiguration
from ....auxiliary.two_dimensional.termination import Termination


class AtomicLayersUniqueRepeatedConfiguration(CrystalLatticePlanesConfiguration, AtomicLayersUniqueRepeatedSchema):
    termination_top: Optional[Termination]
    number_of_repetitions: int
