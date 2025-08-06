from typing import Any, List

from mat3ra.esse.models.materials_category_components.operations.core.combinations.merge import (
    MergeMethodsEnum,
    MergeSchema,
)
from mat3ra.made.tools.build_components.entities.reusable.base_builder import BaseConfigurationPydantic


class MergeConfiguration(MergeSchema, BaseConfigurationPydantic):
    """
    Parameters for merging materials.

    Args:
        merge_components: List of materials or configurations to merge.
        merge_method: Method to use for merging
    """

    merge_components: List[Any]
    merge_method: MergeMethodsEnum
