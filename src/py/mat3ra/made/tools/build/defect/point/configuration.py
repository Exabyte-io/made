from typing import List, Union

from mat3ra.esse.models.materials_category.defective_structures.zero_dimensional.point_defect.base_configuration import (
    PointDefectBaseConfigurationSchema,
)
from mat3ra.esse.models.materials_category.defective_structures.zero_dimensional.point_defect.interstitial import (
    InterstitialPointDefectSchema,
)
from mat3ra.esse.models.materials_category.defective_structures.zero_dimensional.point_defect.substitutional import (
    SubstitutionalPointDefectSchema,
)
from mat3ra.esse.models.materials_category.defective_structures.zero_dimensional.point_defect.vacancy import (
    VacancyPointDefectSchema,
)
from mat3ra.esse.models.materials_category_components.entities.auxiliary.zero_dimensional.point_defect_site import (
    PointDefectSiteSchema,
    AtomSchema,
)
from mat3ra.esse.models.materials_category_components.entities.core.zero_dimensional.vacancy import VacancySchema
from mat3ra.esse.models.materials_category_components.operations.core.combinations.merge import MergeMethodsEnum

from mat3ra.made.material import Material
from mat3ra.made.tools.build.merge.configuration import MergeConfiguration
from mat3ra.made.tools.site import CrystalSite


class PointDefectSite(CrystalSite, PointDefectSiteSchema):
    element: Union[VacancySchema, AtomSchema]


class PointDefectConfiguration(MergeConfiguration, PointDefectBaseConfigurationSchema):
    """
    Configuration for building a point defect by merging materials.

    Args:
        merge_components: List of materials or configurations to merge.
        merge_method: Method to use for merging.
    """

    type: str = "PointDefectConfiguration"
    merge_components: List[Union[Material, PointDefectSite]]


class VacancyDefectConfiguration(PointDefectConfiguration, VacancyPointDefectSchema):
    type: str = "VacancyDefectConfiguration"

    @classmethod
    def from_parameters(cls, crystal: Material, coordinate: List[float], **kwargs):
        point_defect_site = PointDefectSite(
            crystal=crystal,
            element=VacancySchema(chemical_element=VacancySchema.chemical_element.Vac),
            coordinate=coordinate,
        )
        return cls(merge_components=[crystal, point_defect_site], merge_method=MergeMethodsEnum.replace, **kwargs)


class SubstitutionalDefectConfiguration(PointDefectConfiguration, SubstitutionalPointDefectSchema):
    type: str = "SubstitutionalDefectConfiguration"

    @classmethod
    def from_parameters(cls, crystal: Material, coordinate: List[float], element: str, **kwargs):
        substitution_site = PointDefectSite(
            crystal=crystal,
            element=AtomSchema(chemical_element=element),
            coordinate=coordinate,
        )
        return cls(merge_components=[crystal, substitution_site], merge_method=MergeMethodsEnum.replace, **kwargs)


class InterstitialDefectConfiguration(PointDefectConfiguration, InterstitialPointDefectSchema):
    type: str = "InterstitialDefectConfiguration"

    @classmethod
    def from_parameters(cls, crystal: Material, coordinate: List[float], element: str, **kwargs):
        interstitial_site = PointDefectSite(
            crystal=crystal,
            element=AtomSchema(chemical_element=element),
            coordinate=coordinate,
        )
        return cls(merge_components=[crystal, interstitial_site], merge_method=MergeMethodsEnum.add, **kwargs)
