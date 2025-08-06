from mat3ra.esse.models.materials_category.defective_structures.two_dimensional.adatom.configuration import (
    AdatomDefectConfigurationSchema,
    AxisEnum,
)

from .....build_components.entities.reusable.two_dimensional.slab_stack.configuration import SlabStackConfiguration


class AdatomDefectConfiguration(SlabStackConfiguration, AdatomDefectConfigurationSchema):
    """
    Configuration for creating an adatom defect on a slab surface.

    Args:
        stack_components: List containing [slab, isolated_defect, vacuum].
    """

    type: str = "AdatomDefectConfiguration"
    direction: AxisEnum = AxisEnum.z
