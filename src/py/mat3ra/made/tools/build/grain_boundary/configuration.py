from typing import Optional, List

from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum

from ..interface.configuration import InterfaceConfiguration
from ..slab.configurations import SlabConfiguration


class GrainBoundaryConfiguration(InterfaceConfiguration):
    """
    Configuration for creating a grain boundary between two phases.
    """

    direction: AxisEnum = AxisEnum.x

    @property
    def phase_1_configuration(self) -> SlabConfiguration:
        return self.stack_components[0]

    @property
    def phase_2_configuration(self) -> SlabConfiguration:
        return self.stack_components[1]

    @property
    def translation_vector(self) -> List[float]:
        return [0.0, self.xy_shift[0], self.xy_shift[1]]

    @property
    def gap(self) -> float:
        return self.vacuum_configuration.size if len(self.stack_components) > 2 else 0.0

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

        stack_components = [phase_1_configuration, phase_2_configuration]
        return cls(stack_components=stack_components, direction=AxisEnum.x, xy_shift=xy_shift)
