from mat3ra.esse.models.materials_category_components.entities.reusable.two_dimensional.slab_stack_configuration import (  # noqa: E501
    SlabStackConfigurationSchema,
)
from mat3ra.made.material import Material

from .....operations.core.combinations.stack.configuration import StackConfiguration
from ....core.two_dimensional.vacuum.configuration import VacuumConfiguration


class SlabStackConfiguration(StackConfiguration, SlabStackConfigurationSchema):
    """
    Configuration for stacking a slab, an isolated defect and a vacuum layer.

    Args:
        stack_components: List containing [slab, slab_component, vacuum].
        direction: Direction of stacking (default: z).
    """

    type: str = "SlabStackConfiguration"

    @property
    def slab(self) -> Material:
        return self.stack_components[0]

    @property
    def added_component(self) -> Material:
        return self.stack_components[1]

    @property
    def vacuum_configuration(self) -> VacuumConfiguration:
        return self.stack_components[2]
