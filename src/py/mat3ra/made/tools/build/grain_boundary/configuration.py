from typing import Optional, List

from mat3ra.code.array_with_ids import ArrayWithIds
from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
# fmt: off
from mat3ra.esse.models.materials_category.defective_structures.one_dimensional. \
    grain_boundary_linear.configuration import (
    GrainBoundaryLinearConfigurationSchema,
)

from mat3ra.esse.models.materials_category.defective_structures.two_dimensional. \
    grain_boundary_planar.configuration import (
    GrainBoundaryPlanarConfigurationSchema,
)

from ..interface.configuration import InterfaceConfiguration
from ..slab.configurations import SlabConfiguration, SlabStrainedSupercellConfiguration
from ..vacuum.configuration import VacuumConfiguration


# fmt: on


class GrainBoundaryPlanarConfiguration(InterfaceConfiguration, GrainBoundaryPlanarConfigurationSchema):
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
        phase_1_configuration: SlabStrainedSupercellConfiguration,
        phase_2_configuration: SlabStrainedSupercellConfiguration,
        xy_shift: Optional[List[float]] = None,
        gaps: Optional[List[float]] = None,
    ) -> "GrainBoundaryPlanarConfiguration":
        if xy_shift is None:
            xy_shift = [0.0, 0.0]

        stack_components = [phase_1_configuration, phase_2_configuration]
        return cls(
            stack_components=stack_components,
            direction=AxisEnum.z,
            xy_shift=xy_shift,
            gaps=ArrayWithIds.from_values(gaps),
        )


class GrainBoundaryLinearConfiguration(InterfaceConfiguration, GrainBoundaryLinearConfigurationSchema):
    """
    Configuration for creating a linear grain boundary.

    Args:
        stack_components (List): List of configuration objects for grain boundary components.
        direction (AxisEnum): Direction along which to stack components (x or y).
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
