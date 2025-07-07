from typing import Optional, Union, List

from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from mat3ra.esse.models.materials_category.compound_pristine_structures.two_dimensional.interface.configuration import (  # noqa: E501
    InterfaceConfigurationSchema,
)

from mat3ra.made.material import Material
from mat3ra.made.tools.build.vacuum.configuration import VacuumConfiguration
from .. import BaseConfiguration, BaseConfigurationPydantic
from ..slab.configurations import (
    SlabConfiguration,
    SlabStrainedSupercellConfiguration,
    SlabStrainedSupercellWithGapConfiguration,
)


class InterfaceConfiguration(InterfaceConfigurationSchema, BaseConfigurationPydantic):
    # components and their modifiers added in the order they are stacked, from bottom to top
    stack_components: List[
        Union[
            SlabStrainedSupercellConfiguration,
            SlabStrainedSupercellWithGapConfiguration,
            VacuumConfiguration,
        ]
    ]
    direction: AxisEnum = AxisEnum.z
    xy_shift: List[float] = InterfaceConfigurationSchema.model_fields["xy_shift"].default  # in Angstroms

    @property
    def substrate_configuration(self) -> SlabConfiguration:
        return self.stack_components[0]

    @property
    def film_configuration(self) -> SlabConfiguration:
        return self.stack_components[1]

    @property
    def vacuum_configuration(self) -> VacuumConfiguration:
        if len(self.stack_components) > 2:
            return self.stack_components[2]
        return VacuumConfiguration(
            size=0.0, crystal=self.film_configuration.atomic_layers.crystal, direction=self.direction
        )


class TwistedNanoribbonsInterfaceConfiguration(BaseConfiguration):
    """
    Configuration for creating a twisted interface between two nanoribbons with specified twist angle.

    Args:
        stack_components (List[Material]): List of two nanoribbon materials.
        angle (float): Twist angle in degrees for provenance.
    """

    stack_components: List[Material]
    angle: float = 0.0

    @property
    def nanoribbon1(self) -> Material:
        return self.stack_components[0]

    @property
    def nanoribbon2(self) -> Material:
        return self.stack_components[1]
