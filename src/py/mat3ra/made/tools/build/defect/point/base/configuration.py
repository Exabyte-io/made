from typing import Union, List

from mat3ra.esse.models.materials_category.defective_structures.zero_dimensional.point_defect.base_configuration import (
    PointDefectBaseConfigurationSchema,
)

from mat3ra.made.material import Material
from mat3ra.made.tools.build.defect.point.defect_site.configuration import PointDefectSiteConfiguration
from mat3ra.made.tools.build.merge.configuration import MergeConfiguration


class PointDefectConfiguration(MergeConfiguration, PointDefectBaseConfigurationSchema):
    """
    Configuration for building a point defect by merging materials.

    Args:
        merge_components: List of materials or configurations to merge.
        merge_method: Method to use for merging.
    """

    type: str = "PointDefectConfiguration"
    merge_components: List[Union[Material, PointDefectSiteConfiguration]]
