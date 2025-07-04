from typing import List

from mat3ra.esse.models.materials_category.defective_structures.zero_dimensional.point_defect.base_configuration import (
    PointDefectBaseConfigurationSchema,
)
from mat3ra.esse.models.materials_category.defective_structures.zero_dimensional.point_defect.substitutional import (
    SubstitutionalPointDefectSchema,
    PointDefectSiteSchema,
)
from mat3ra.esse.models.materials_category.defective_structures.zero_dimensional.point_defect.vacancy import (
    VacancyPointDefectSchema,
)
from mat3ra.esse.models.materials_category_components.entities.core.zero_dimensional.vacancy import VacancySchema

from mat3ra.made.material import Material
from mat3ra.made.tools.build.merge.configuration import MergeConfiguration
from mat3ra.made.tools.site import CrystalSite


class PointDefectSite(PointDefectSiteSchema, CrystalSite):
    pass


class PointDefectConfiguration(MergeConfiguration, PointDefectBaseConfigurationSchema):
    """
    Configuration for building a point defect by merging materials.

    Args:
        merge_components: List of materials or configurations to merge.
        merge_method: Method to use for merging.
    """

    type: str = "PointDefectConfiguration"
    merge_components: List[Material, PointDefectSite]


class VacancyDefectConfiguration(VacancyPointDefectSchema, PointDefectConfiguration):
    """
    Configuration for building a vacancy defect by merging materials.

    Args:
        merge_components: List of materials or configurations to merge.
        merge_method: Method to use for merging.
    """

    type: str = "VacancyDefectConfiguration"
    merge_components: List[Material, VacancySchema]


class SubstitutionalDefectConfiguration(SubstitutionalPointDefectSchema, PointDefectConfiguration):
    """
    Configuration for building a substitutional defect by merging materials.

    Args:
        merge_components: List of materials or configurations to merge.
        merge_method: Method to use for merging.
    """

    type: str = "SubstitutionalDefectConfiguration"
    merge_components: List[Material, PointDefectSite]


class InterstitialDefectConfiguration(PointDefectConfiguration):
    """
    Configuration for building an interstitial defect by merging materials.

    Args:
        merge_components: List of materials or configurations to merge.
        merge_method: Method to use for merging.
    """

    type: str = "InterstitialDefectConfiguration"
    merge_components: List[Material, PointDefectSite]
