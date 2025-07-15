from typing import List, Union

from mat3ra.esse.models.materials_category.defective_structures.zero_dimensional.point_defect.base_configuration import (
    PointDefectBaseConfigurationSchema,
)
from mat3ra.esse.models.materials_category_components.operations.core.combinations.merge import MergeMethodsEnum

from mat3ra.made.material import Material
from mat3ra.made.tools.build.defect.point.configuration import PointDefectConfiguration
from mat3ra.made.tools.build.merge.configuration import MergeConfiguration


class PairDefectConfiguration(MergeConfiguration, PointDefectBaseConfigurationSchema):
    """
    Configuration for building a pair defect by merging materials.

    Args:
        merge_components: List containing [crystal, primary_defect_configuration, secondary_defect_configuration].
        merge_method: Method to use for merging.
    """

    type: str = "PairDefectConfiguration"
    merge_components: List[Union[Material, PointDefectConfiguration]]
    merge_method: MergeMethodsEnum = MergeMethodsEnum.REPLACE

    @classmethod
    def from_parameters(
        cls,
        crystal: Material,
        primary_defect_configuration: PointDefectConfiguration,
        secondary_defect_configuration: PointDefectConfiguration,
        **kwargs,
    ):
        return cls(merge_components=[crystal, primary_defect_configuration, secondary_defect_configuration], **kwargs)
