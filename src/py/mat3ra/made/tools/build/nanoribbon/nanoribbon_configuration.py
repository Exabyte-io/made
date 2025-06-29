from typing import List, Union
from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from mat3ra.made.tools.build import BaseConfigurationPydantic
from ..stack.configuration import StackConfiguration
from ..vacuum.configuration import VacuumConfiguration
from .nano_tape_configuration import NanoTapeConfiguration


class NanoribbonConfiguration(StackConfiguration):
    """
    Configuration for building a nanoribbon from a nanotape.
    Nanoribbon = [NanoTape, vacuum] stacked on X direction.

    Args:
        stack_components: List of configuration objects for nanoribbon components.
        direction: Direction along which to stack components.
    """

    type: str = "NanoribbonConfiguration"
    stack_components: List[Union[NanoTapeConfiguration, VacuumConfiguration]]
    direction: AxisEnum = AxisEnum.x
    use_rectangular_lattice: bool = True

    @property
    def nanotape(self):
        """Get the nanotape configuration component."""
        return self.stack_components[0]

    @property
    def vacuum_configuration(self) -> VacuumConfiguration:
        """Get the vacuum configuration component."""
        return self.stack_components[1] 