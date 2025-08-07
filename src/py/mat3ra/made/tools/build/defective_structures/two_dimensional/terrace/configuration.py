from typing import List

from mat3ra.esse.models.materials_category.defective_structures.two_dimensional.terrace.configuration import (
    TerraceDefectConfigurationSchema,
)

from .....build_components.entities.reusable.two_dimensional.slab_stack.configuration import SlabStackConfiguration


class TerraceDefectConfiguration(SlabStackConfiguration, TerraceDefectConfigurationSchema):
    cut_direction: List[int] = [1, 0, 0]
