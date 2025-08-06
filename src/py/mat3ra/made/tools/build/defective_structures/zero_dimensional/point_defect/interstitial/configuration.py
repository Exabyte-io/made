from typing import List, Union

from mat3ra.esse.models.materials_category.defective_structures.zero_dimensional.point_defect.interstitial import (
    InterstitialPointDefectSchema,
)
from mat3ra.esse.models.materials_category_components.entities.core.zero_dimensional.atom import AtomSchema
from mat3ra.esse.models.materials_category_components.operations.core.combinations.merge import MergeMethodsEnum

from mat3ra.made.material import Material
from ......build_components import MaterialWithBuildMetadata
from ..base.configuration import PointDefectConfiguration
from ......build_components.entities.auxiliary.zero_dimensional.point_defect_site.configuration import (
    PointDefectSiteConfiguration,
)


class InterstitialDefectConfiguration(PointDefectConfiguration, InterstitialPointDefectSchema):
    type: str = "InterstitialDefectConfiguration"

    @classmethod
    def from_parameters(
        cls, crystal: Union[Material, MaterialWithBuildMetadata], coordinate: List[float], element: str, **kwargs
    ):
        interstitial_site = PointDefectSiteConfiguration(
            crystal=crystal,
            element=AtomSchema(chemical_element=element),
            coordinate=coordinate,
        )
        return cls(merge_components=[crystal, interstitial_site], merge_method=MergeMethodsEnum.ADD, **kwargs)
