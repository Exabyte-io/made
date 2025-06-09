from typing import List, Union

from mat3ra.code.entity import InMemoryEntityPydantic
from mat3ra.esse.models.materials_category_components.entities.auxiliary.two_dimensional.miller_indices import (
    MillerIndicesSchema,
)
from mat3ra.esse.models.materials_category_components.entities.reusable.two_dimensional.crystal_lattice_planes import (
    CrystalLatticePlanesSchema,
)
from pydantic import BaseModel

from mat3ra.made.material import Material
from .. import BaseConfigurationPydantic

from ..stack.configuration import StackConfiguration
from ..vacuum.configuration import VacuumConfiguration
from ...analyze.termination import Termination


class CrystalLatticePlanesConfiguration(CrystalLatticePlanesSchema, BaseConfigurationPydantic):
    crystal: Material


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
