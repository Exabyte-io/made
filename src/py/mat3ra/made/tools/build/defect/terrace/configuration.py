from mat3ra.esse.models.materials_category.defective_structures.two_dimensional.terrace.configuration import (
    TerraceDefectConfigurationSchema,
)
from mat3ra.made.tools.build.defect.slab.configuration import SlabStackConfiguration


class TerraceDefectConfiguration(SlabStackConfiguration, TerraceDefectConfigurationSchema):
    cut_direction: list[int] = [1, 0, 0]
