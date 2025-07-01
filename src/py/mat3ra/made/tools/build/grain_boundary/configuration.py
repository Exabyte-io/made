from typing import Optional, List, Union

from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum

from .. import BaseConfiguration
from ..slab.configurations import SlabConfiguration
from ..slab.entities import Termination
from ..interface.configuration import TwistedInterfaceConfiguration
from ..stack.configuration import StackConfiguration
from ..vacuum.configuration import VacuumConfiguration


class SlabGrainBoundaryConfiguration(BaseConfiguration):
    """
    Configuration for a grain boundary between two phases with different surfaces facing each other.

    Attributes:
        phase_1_configuration (SlabConfiguration): The configuration of the first phase.
        phase_2_configuration (SlabConfiguration): The configuration of the second phase.
        phase_1_termination (Termination): The termination of the first phase.
        phase_2_termination (Termination): The termination of the second phase.
        gap (float): The gap between the two phases, in Angstroms.
        slab_configuration (SlabConfiguration): The configuration of the grain boundary slab.
        slab_termination (Optional[Termination]): The termination of the grain boundary slab.
    """

    phase_1_configuration: SlabConfiguration
    phase_2_configuration: SlabConfiguration
    phase_1_termination: Termination
    phase_2_termination: Termination
    gap: float = 3.0
    slab_configuration: SlabConfiguration
    slab_termination: Optional[Termination] = None

    @property
    def _json(self):
        return {
            "type": self.__class__.__name__,
            "phase_1_configuration": self.phase_1_configuration.to_json(),
            "phase_2_configuration": self.phase_2_configuration.to_json(),
            "phase_1_termination": str(self.phase_1_termination),
            "phase_2_termination": str(self.phase_2_termination),
            "gap": self.gap,
            "slab_configuration": self.slab_configuration.to_json(),
        }


class GrainBoundaryLinearConfiguration(StackConfiguration):
    """
    Configuration for creating a linear grain boundary.

    Args:
        stack_components (List): List of configuration objects for grain boundary components.
        direction (AxisEnum): Direction along which to stack components (x or y).
        gap (float): The gap between the two phases.
    """

    type: str = "GrainBoundaryLinearConfiguration"
    stack_components: List[Union[SlabConfiguration, VacuumConfiguration]]
    direction: AxisEnum = AxisEnum.x
    gap: float = 0.0

    @property
    def _json(self):
        return {
            "type": self.get_cls_name(),
            "stack_components": [comp.to_json() for comp in self.stack_components],
            "direction": self.direction.value,
            "gap": self.gap,
        }
