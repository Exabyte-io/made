from typing import List, Optional

from mat3ra.code.array_with_ids import ArrayWithIds
from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
# fmt: off
from mat3ra.esse.models.materials_category.defective_structures.two_dimensional. \
    grain_boundary_planar.configuration import (
    GrainBoundaryPlanarConfigurationSchema,
)
# fmt: on

from mat3ra.made.tools.build.interface import InterfaceConfiguration
from mat3ra.made.tools.build.slab.slab.configuration import SlabConfiguration
from mat3ra.made.tools.build.slab.strained_supercell_slab.configuration import SlabStrainedSupercellConfiguration
from mat3ra.made.tools.build.vacuum.configuration import VacuumConfiguration


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
