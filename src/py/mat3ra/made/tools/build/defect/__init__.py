from typing import Union, Optional

from mat3ra.made.material import Material
from .builders import VacancyBuilder, SubstitutionBuilder, InterstitialBuilder, PointDefectBuilderParameters
from .configuration import VacancyConfiguration, SubstitutionConfiguration, InterstitialConfiguration
from .enums import PointDefectTypeEnum


def DEFECT_BUILDER_MAP(builder_parameters):
    return {
        PointDefectTypeEnum.VACANCY: VacancyBuilder(builder_parameters),
        PointDefectTypeEnum.SUBSTITUTION: SubstitutionBuilder(builder_parameters),
        PointDefectTypeEnum.INTERSTITIAL: InterstitialBuilder(builder_parameters),
    }


def create_defect(
    material: Material,
    configuration: Union[VacancyConfiguration, SubstitutionConfiguration, InterstitialConfiguration],
    builder_parameters: Optional[PointDefectBuilderParameters] = None,
) -> Material:
    """
    Return a material with a selected defect added.

    Args:
        material: The material to which the defect will be added.
        configuration: The configuration of the defect to be added.

    Returns:
        The material with the defect added.
    """
    builder = DEFECT_BUILDER_MAP(builder_parameters)[configuration.defect_type]

    configuration.crystal = material
    return builder.get_material(configuration) if builder else material
