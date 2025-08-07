from typing import Union, List

from mat3ra.esse.models.materials_category.defective_structures.zero_dimensional.point_defect.vacancy import (
    VacancyPointDefectSchema,
)
from mat3ra.esse.models.materials_category_components.entities.core.zero_dimensional.vacancy import VacancySchema
from mat3ra.esse.models.materials_category_components.operations.core.combinations.merge import MergeMethodsEnum

from mat3ra.made.material import Material
from ..base.configuration import PointDefectConfiguration
from ......build_components import MaterialWithBuildMetadata
from ......build_components.entities.auxiliary.zero_dimensional.point_defect_site.configuration import (
    PointDefectSiteConfiguration,
)


class VacancyDefectConfiguration(PointDefectConfiguration, VacancyPointDefectSchema):
    type: str = "VacancyDefectConfiguration"

    @classmethod
    def from_parameters(cls, crystal: Union[Material, MaterialWithBuildMetadata], coordinate: List[float], **kwargs):
        point_defect_site = PointDefectSiteConfiguration(
            crystal=crystal,
            element=VacancySchema(),
            coordinate=coordinate,
        )
        return cls(merge_components=[crystal, point_defect_site], merge_method=MergeMethodsEnum.REPLACE, **kwargs)
