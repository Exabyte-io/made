from typing import List, Union

from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from mat3ra.esse.models.materials_category.pristine_structures.two_dimensional.nanoribbon import (
    NanoribbonConfigurationSchema,
)

from ..nanotape.configuration import NanoTapeConfiguration
from ..stack.configuration import StackConfiguration
from ..vacuum.configuration import VacuumConfiguration


class NanoribbonConfiguration(NanoribbonConfigurationSchema, StackConfiguration):
    """
    Configuration for building a nanoribbon from a nanotape.
    Nanoribbon = [NanoTape, vacuum] stacked in X or Y direction.

    Args:
        stack_components: List of configuration objects for nanoribbon components.
        direction: Direction along which to stack components.
    """

    type: str = "NanoribbonConfiguration"
    stack_components: List[Union[NanoTapeConfiguration, VacuumConfiguration]]
    direction: AxisEnum = AxisEnum.x

    @property
    def nanotape(self):
        return self.stack_components[0]

    @property
    def vacuum_configuration(self) -> VacuumConfiguration:
        return self.stack_components[1]
