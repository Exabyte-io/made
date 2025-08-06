from typing import List, Union

from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from mat3ra.esse.models.materials_category.pristine_structures.two_dimensional.nanotape import (
    NanoTapeConfigurationSchema,
)

from mat3ra.made.tools.build_components.entities.core.two_dimensional.vacuum.configuration import VacuumConfiguration
from mat3ra.made.tools.build_components.operations.core.combinations.stack.configuration import StackConfiguration
from mat3ra.made.tools.build_components.operations.core.modifications.repeat import \
    CrystalLatticeLinesUniqueRepeatedConfiguration


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
