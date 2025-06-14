from typing import List, Union

from mat3ra.esse.models.materials_category_components.entities.reusable.two_dimensional.crystal_lattice_planes import (
    CrystalLatticePlanesSchema,
)

from mat3ra.made.material import Material
from mat3ra.made.tools.build.slab.entities import Termination
from .. import BaseConfigurationPydantic
from ..stack.configuration import StackConfiguration
from ..vacuum.configuration import VacuumConfiguration


class CrystalLatticePlanesConfiguration(CrystalLatticePlanesSchema, BaseConfigurationPydantic):
    crystal: Material


class AtomicLayersUnique(CrystalLatticePlanesConfiguration):
    pass


class AtomicLayersUniqueRepeatedConfiguration(AtomicLayersUnique):
    termination_top: Termination
    number_of_repetitions: int = 1


class SlabConfiguration(StackConfiguration):
    type: str = "SlabConfiguration"
    stack_components: List[
        Union[AtomicLayersUnique, AtomicLayersUniqueRepeatedConfiguration, VacuumConfiguration]  # No Materials!
    ]

    @property
    def atomic_layers(self):
        return self.stack_components[0]

    @property
    def vacuum_configuration(self) -> VacuumConfiguration:
        return self.stack_components[1]
