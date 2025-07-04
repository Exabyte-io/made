from typing import List, Union

from mat3ra.esse.models.element import ElementSchema
from mat3ra.esse.models.materials_category.defective_structures.zero_dimensional.point_defect.base_configuration import (
    PointDefectBaseConfigurationSchema,
)
from mat3ra.esse.models.materials_category.defective_structures.zero_dimensional.point_defect.substitutional import (
    SubstitutionalPointDefectSchema,
)
from mat3ra.esse.models.materials_category.defective_structures.zero_dimensional.point_defect.vacancy import (
    VacancyPointDefectSchema,
)
from mat3ra.esse.models.materials_category_components.entities.auxiliary.zero_dimensional.point_defect_site import (
    PointDefectSiteSchema,
)
from mat3ra.esse.models.materials_category_components.entities.core.zero_dimensional.vacancy import VacancySchema
from mat3ra.esse.models.materials_category_components.operations.core.combinations.merge import MergeMethodsEnum

from mat3ra.made.material import Material
from mat3ra.made.tools.build.merge.configuration import MergeConfiguration
from mat3ra.made.tools.site import CrystalSite


class PointDefectSite(CrystalSite, PointDefectSiteSchema):
    element: Union[VacancySchema, ElementSchema]


class PointDefectConfiguration(MergeConfiguration, PointDefectBaseConfigurationSchema):
    """
    Configuration for building a point defect by merging materials.

    Args:
        merge_components: List of materials or configurations to merge.
        merge_method: Method to use for merging.
    """

    type: str = "PointDefectConfiguration"
    merge_components: List[Union[Material, PointDefectSite]]


class VacancyDefectConfiguration(VacancyPointDefectSchema, PointDefectConfiguration):
    """
    Configuration for building a vacancy defect by merging materials.

    Args:
        merge_components: List of materials or configurations to merge.
        merge_method: Method to use for merging.
    """

    type: str = "VacancyDefectConfiguration"

    @classmethod
    def from_parameters(cls, host_material: Material, coordinate: List[float], **kwargs):

        point_defect_site = PointDefectSite(
            element=VacancySchema(chemical_element=VacancySchema.chemical_element.Vac),
            coordinate=coordinate,
        )

        return cls(merge_components=[host_material, point_defect_site], merge_method=MergeMethodsEnum.replace, **kwargs)


class SubstitutionalDefectConfiguration(SubstitutionalPointDefectSchema, PointDefectConfiguration):
    """
    Configuration for building a substitutional defect by merging materials.

    Args:
        merge_components: List of materials or configurations to merge.
        merge_method: Method to use for merging.
    """

    type: str = "SubstitutionalDefectConfiguration"
    merge_components: List[Union[Material, PointDefectSite]]

    @classmethod
    def from_parameters(cls, host_material: Material, coordinate: List[float], element: str, **kwargs):
        substitution_site = PointDefectSite(element=element, coordinate=coordinate)

        return cls(merge_components=[host_material, substitution_site], merge_method=MergeMethodsEnum.replace, **kwargs)


class InterstitialDefectConfiguration(PointDefectConfiguration):
    """
    Configuration for building an interstitial defect by merging materials.

    Args:
        merge_components: List of materials or configurations to merge.
        merge_method: Method to use for merging.
    """

    type: str = "InterstitialDefectConfiguration"
    merge_components: List[Union[Material, PointDefectSite]]

    @classmethod
    def from_parameters(cls, host_material: Material, coordinate: List[float], element: str, **kwargs):

        # TODO: convert str to correct ElementSchema
        interstitial_site = PointDefectSite(
            element=element,
            coordinate=coordinate,
        )

        return cls(merge_components=[host_material, interstitial_site], merge_method=MergeMethodsEnum.add, **kwargs)
