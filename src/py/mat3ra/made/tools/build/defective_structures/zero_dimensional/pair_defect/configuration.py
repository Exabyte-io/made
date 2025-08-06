from typing import List, Union

# fmt: off
from mat3ra.esse.models.materials_category.defective_structures.zero_dimensional. \
    point_defect.base_configuration import (
    PointDefectBaseConfigurationSchema,
)
# fmt: on
from mat3ra.esse.models.materials_category_components.operations.core.combinations.merge import MergeMethodsEnum

from mat3ra.made.material import Material
from ..point_defect.base.configuration import PointDefectConfiguration
from .....build_components import MaterialWithBuildMetadata
from .....build_components.entities.auxiliary.zero_dimensional.point_defect_site.configuration import (
    PointDefectSiteConfiguration,
)
from .....build_components.operations.core.combinations.merge.configuration import MergeConfiguration


class PairDefectConfiguration(MergeConfiguration, PointDefectBaseConfigurationSchema):
    """
    Configuration for building a pair defect by merging materials.

    Args:
        merge_components: List containing [crystal, primary_defect_configuration, secondary_defect_configuration].
        merge_method: Method to use for merging.
    """

    type: str = "PairDefectConfiguration"
    merge_components: List[Union[Material, PointDefectSiteConfiguration]]
    merge_method: MergeMethodsEnum = MergeMethodsEnum.REPLACE

    @classmethod
    def from_parameters(
        cls,
        crystal: Union[Material, MaterialWithBuildMetadata],
        primary_defect_configuration: PointDefectConfiguration,
        secondary_defect_configuration: PointDefectConfiguration,
        **kwargs,
    ):
        primary_defect_site_configuration = primary_defect_configuration.merge_components[1]
        secondary_defect_site_configuration = secondary_defect_configuration.merge_components[1]
        return cls(
            merge_components=[crystal, primary_defect_site_configuration, secondary_defect_site_configuration], **kwargs
        )
