from typing import List
from mat3ra.esse.models.materials_category.defective_structures.two_dimensional.terrace.configuration import (
    TerraceDefectConfigurationSchema,
)
from ..slab.configuration import SlabStackConfiguration


class TerraceDefectConfiguration(SlabStackConfiguration, TerraceDefectConfigurationSchema):
    cut_direction: List[int] = [1, 0, 0]
