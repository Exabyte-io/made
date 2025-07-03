from typing import Optional, List

from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from mat3ra.esse.models.materials_category.defective_structures.one_dimensional.grain_boundary_linear.configuration import (
    GrainBoundaryLinearConfigurationSchema,
)
# fmt: off
from mat3ra.esse.models.materials_category.defective_structures.two_dimensional. \
    grain_boundary_planar.configuration import (
    GrainBoundaryPlanarConfigurationSchema as GrainBoundarySchema,
)

from ..interface.configuration import InterfaceConfiguration
from ..slab.configurations import SlabConfiguration, SlabStrainedSupercellWithGapConfiguration
from ..vacuum.configuration import VacuumConfiguration


# fmt: on


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


class GrainBoundaryLinearConfiguration(InterfaceConfiguration, GrainBoundaryLinearConfigurationSchema):
    """
    Configuration for creating a linear grain boundary.

    Args:
        stack_components (List): List of configuration objects for grain boundary components.
        direction (AxisEnum): Direction along which to pypstack components (x or y).
    """

    type: str = "GrainBoundaryLinearConfiguration"
    direction: AxisEnum = AxisEnum.x
    actual_angle: Optional[float] = None

    @property
    def phase_1_configuration(self) -> SlabConfiguration:
        return self.stack_components[0]

    @property
    def phase_2_configuration(self) -> SlabConfiguration:
        return self.stack_components[1]

    @property
    def substrate_configuration(self) -> SlabConfiguration:
        return self.phase_1_configuration

    @property
    def film_configuration(self) -> SlabConfiguration:
        return self.phase_2_configuration
