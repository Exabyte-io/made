from typing import Optional, Union

from mat3ra.utils.factory import BaseFactory
from mat3ra.made.material import Material

from .builders import (
    PointDefectBuilderParameters,
    SlabDefectBuilderParameters,
    AdatomSlabDefectBuilder,
    EquidistantAdatomSlabDefectBuilder,
)
from .configuration import PointDefectConfiguration, AdatomSlabDefectConfiguration
from .enums import PointDefectTypeEnum, SlabDefectTypeEnum


class DefectBuilderFactory(BaseFactory):
    __class_registry__ = {
        PointDefectTypeEnum.VACANCY: "mat3ra.made.tools.build.defect.builders.VacancyPointDefectBuilder",
        PointDefectTypeEnum.SUBSTITUTION: "mat3ra.made.tools.build.defect.builders.SubstitutionPointDefectBuilder",
        PointDefectTypeEnum.INTERSTITIAL: "mat3ra.made.tools.build.defect.builders.InterstitialPointDefectBuilder",
    }


def create_defect(
    configuration: Union[PointDefectConfiguration, AdatomSlabDefectConfiguration],
    builder_parameters: Union[PointDefectBuilderParameters, SlabDefectBuilderParameters, None] = None,
) -> Material:
    """
    Return a material with a selected defect added.

    Args:
        configuration: The configuration of the defect to be added.
        builder_parameters: The parameters to be used by the defect builder.

    Returns:
        The material with the defect added.
    """
    BuilderClass = DefectBuilderFactory.get_class_by_name(configuration.defect_type)
    builder = BuilderClass(builder_parameters)

    return builder.get_material(configuration) if builder else configuration.crystal


def create_slab_defect(
    configuration: Union[AdatomSlabDefectConfiguration],
    builder: Optional[Union[AdatomSlabDefectBuilder, EquidistantAdatomSlabDefectBuilder]] = None,
) -> Material:
    """
    Return a material with a selected slab defect added.

    Args:
        configuration: The configuration of the defect to be added.
        builder: The builder to be used to create the defect.

    Returns:
        The material with the defect added.
    """
    if builder is None:
        builder = AdatomSlabDefectBuilder(build_parameters=SlabDefectBuilderParameters())
    return builder.get_material(configuration)
