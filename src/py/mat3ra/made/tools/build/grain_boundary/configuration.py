from typing import Optional, List

from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
# fmt: off
from mat3ra.esse.models.materials_category.defective_structures.two_dimensional.\
    grain_boundary_planar.configuration import (
    GrainBoundaryPlanarConfigurationSchema as GrainBoundarySchema,
)
# fmt: on

from ..interface.configuration import InterfaceConfiguration
from ..slab.configurations import SlabConfiguration, SlabStrainedSupercellWithGapConfiguration
from ..vacuum.configuration import VacuumConfiguration


class GrainBoundaryPlanarConfiguration(InterfaceConfiguration, GrainBoundarySchema):
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
    ) -> "GrainBoundaryPlanarConfiguration":
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
