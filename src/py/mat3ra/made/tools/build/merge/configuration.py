from typing import List, Any

from mat3ra.esse.models.materials_category_components.operations.core.combinations.merge import (
    MergeSchema,
    MergeMethodsEnum,
)

from mat3ra.made.tools.build import BaseConfigurationPydantic


class MergeConfiguration(MergeSchema, BaseConfigurationPydantic):
    """
    Parameters for merging materials.

    Args:
        merge_components: List of materials or configurations to merge.
        merge_method: Method to use for merging
    """

    merge_components: List[Any]
    merge_method: MergeMethodsEnum
