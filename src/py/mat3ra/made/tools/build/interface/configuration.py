from typing import Union, List, Optional

from mat3ra.esse.models.materials_category.compound_pristine_structures.two_dimensional.interface.configuration import (  # noqa: E501
    InterfaceConfigurationSchema,
)

from mat3ra.made.tools.build.vacuum.configuration import VacuumConfiguration
from ..slab.configurations import (
    SlabConfiguration,
    SlabStrainedSupercellConfiguration,
    SlabStrainedSupercellWithGapConfiguration,
)
from ..stack.configuration import StackConfiguration


class InterfaceConfiguration(StackConfiguration, InterfaceConfigurationSchema):
    # components and their modifiers added in the order they are stacked, from bottom to top
    stack_components: List[
        Union[
            SlabStrainedSupercellConfiguration,
            SlabStrainedSupercellWithGapConfiguration,
            VacuumConfiguration,
        ]
    ]
    xy_shift: List[float] = InterfaceConfigurationSchema.model_fields["xy_shift"].default  # in Angstroms

    @property
    def substrate_configuration(self) -> SlabConfiguration:
        return self.stack_components[0]

    @property
    def film_configuration(self) -> SlabConfiguration:
        return self.stack_components[1]

    @property
    def vacuum_configuration(self) -> Optional[VacuumConfiguration]:
        if len(self.stack_components) > 2:
            return self.stack_components[2]
        return None


class TwistedNanoribbonsInterfaceConfiguration(InterfaceConfiguration):
    """
    Configuration for creating a twisted interface between two nanoribbons with specified twist angle.

    Args:
        stack_components (List[SlabConfiguration]): List of two nanoribbons as slab configurations.
        angle (float): Twist angle in degrees for provenance.
    """

    angle: float = 0.0
