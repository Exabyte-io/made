from typing import Optional, Union, List

from mat3ra.made.material import Material

from .builders import (
    PointDefectBuilderParameters,
    SlabDefectBuilderParameters,
    AdatomSlabDefectBuilder,
    EquidistantAdatomSlabDefectBuilder,
    CrystalSiteAdatomSlabDefectBuilder,
    IslandSlabDefectBuilder,
    TerraceSlabDefectBuilder,
)
from .configuration import PointDefectConfiguration, AdatomSlabPointDefectConfiguration, IslandSlabDefectConfiguration
from .enums import PointDefectTypeEnum
from .factories import DefectBuilderFactory


def create_defect(
    configuration: Union[PointDefectConfiguration, AdatomSlabPointDefectConfiguration],
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


def create_defects(
    configurations: Union[List[PointDefectConfiguration], List[AdatomSlabPointDefectConfiguration]],
    builder_parameters: Union[PointDefectBuilderParameters, SlabDefectBuilderParameters, None] = None,
) -> Material:
    """
    Return a material with accumulated defects added.

    Args:
        configurations: The list of configurations of the defect to be added. The defects will be added in that order.
        builder_parameters: The parameters to be used by the defect builder.

    Returns:
        The material with the defects added.
    """
    material_with_defect = None

    for configuration in configurations:
        if material_with_defect:
            configuration.crystal = material_with_defect

        BuilderClass = DefectBuilderFactory.get_class_by_name(configuration.defect_type)
        builder = BuilderClass(builder_parameters)
        material_with_defect = builder.get_material(configuration) if builder else configuration.crystal

    return material_with_defect


def create_slab_defect(
    configuration: Union[AdatomSlabPointDefectConfiguration, IslandSlabDefectConfiguration, TerraceSlabDefectBuilder],
    builder: Optional[
        Union[
            AdatomSlabDefectBuilder,
            EquidistantAdatomSlabDefectBuilder,
            CrystalSiteAdatomSlabDefectBuilder,
            IslandSlabDefectBuilder,
            TerraceSlabDefectBuilder,
        ]
    ] = None,
) -> Material:
    """
    Return a material with a selected slab defect added.

    Args:
        configuration: The configuration of the defect to be added.
        builder: Optional builder to be used to create the defect.
            If None, builder will be selected based on configuration.

    Returns:
        The material with the defect added.
    """
    if builder is None:
        if configuration.defect_type == PointDefectTypeEnum.ADATOM:
            builder_key = f"adatom:{configuration.placement_method.lower()}"
        else:
            builder_key = configuration.defect_type.lower()

        BuilderClass = DefectBuilderFactory.get_class_by_name(builder_key)
        builder = BuilderClass(build_parameters=SlabDefectBuilderParameters())

    return builder.get_material(configuration)
