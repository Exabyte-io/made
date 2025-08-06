from typing import Any, List

from mat3ra.esse.models.materials_category_components.operations.core.combinations.merge import (
    MergeMethodsEnum,
    MergeSchema,
)
from mat3ra.made.tools.build_components.entities.reusable.three_dimensional.crystal_lattice_base.base_configuration_pydantic import (
    BaseConfigurationPydantic,
)


class MergeConfiguration(MergeSchema, BaseConfigurationPydantic):
    """
    Parameters for merging materials.

    Args:
        merge_components: List of materials or configurations to merge.
        merge_method: Method to use for merging
    """

    merge_components: List[Any]
    merge_method: MergeMethodsEnum
