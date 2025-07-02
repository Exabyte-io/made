from typing import Optional, List, Union

from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum

from .. import BaseConfiguration
from ..slab.configurations import SlabConfiguration, SlabStrainedSupercellWithGapConfiguration
from ..slab.entities import Termination
from ..stack.configuration import StackConfiguration
from ..vacuum.configuration import VacuumConfiguration

from ..interface.configuration import InterfaceConfiguration


class GrainBoundaryConfiguration(InterfaceConfiguration):
    """
    Configuration for creating a grain boundary between two phases.
    """

    direction: AxisEnum = AxisEnum.z

    @property
    def phase_1_configuration(self) -> SlabConfiguration:
        return self.stack_components[0]

    @property
    def phase_2_configuration(self) -> SlabConfiguration:
        return self.stack_components[1]

    @property
    def vacuum_configuration(self) -> VacuumConfiguration:
        if len(self.stack_components) > 2:
            return self.stack_components[2]
        return VacuumConfiguration(
            size=0.0, crystal=self.film_configuration.atomic_layers.crystal, direction=self.direction
        )

    @classmethod
    def from_parameters(
        cls,
        phase_1_configuration: SlabConfiguration,
        phase_2_configuration: SlabConfiguration,
        xy_shift: Optional[List[float]] = None,
        gap: Optional[float] = None,
    ) -> "GrainBoundaryConfiguration":
        if xy_shift is None:
            xy_shift = [0.0, 0.0]

        if gap and gap > 0:
            phase_1_config = SlabStrainedSupercellWithGapConfiguration(**phase_1_configuration.to_dict(), gap=gap)
            phase_2_config = SlabStrainedSupercellWithGapConfiguration(**phase_2_configuration.to_dict(), gap=gap)
        else:
            phase_1_config = phase_1_configuration
            phase_2_config = phase_2_configuration

        stack_components = [phase_1_config, phase_2_config]
        return cls(stack_components=stack_components, direction=AxisEnum.z, xy_shift=xy_shift)

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
    gap: float = 3.0
