from typing import List, Union

from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from mat3ra.esse.models.materials_category.pristine_structures.two_dimensional.nanotape import (
    NanoTapeConfigurationSchema,
)

from ..lattice_lines.configuration import (
    CrystalLatticeLinesUniqueRepeatedConfiguration,
)
from ..stack.configuration import StackConfiguration
from ..vacuum.configuration import VacuumConfiguration


class NanoTapeConfiguration(NanoTapeConfigurationSchema, StackConfiguration):
    """
    Configuration for building a nanotape from crystal lattice lines.
    NanoTape = [CrystalLatticeLinesUniqueRepeatedConfiguration, vacuum] stacked in X or Y direction.

    Args:
        stack_components: List of configuration objects for nanotape components.
        direction: Direction along which to stack components.
    """

    type: str = "NanoTapeConfiguration"
    stack_components: List[Union[CrystalLatticeLinesUniqueRepeatedConfiguration, VacuumConfiguration]]
    direction: AxisEnum = AxisEnum.y

    @property
    def lattice_lines(self):
        return self.stack_components[0]

    @property
    def vacuum_configuration(self) -> VacuumConfiguration:
        return self.stack_components[1]
