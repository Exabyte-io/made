from typing import Optional, List

from mat3ra.made.material import Material
from mat3ra.made.tools.build.defect.enums import PointDefectTypeEnum, AtomPlacementMethodEnum
from mat3ra.made.tools.build.merge.builders import PointDefectBuilder, MergeBuilder, MergeBuilderParameters

__all__ = [
    "MergeBuilder",
    "MergeBuilderParameters",
    "PointDefectBuilder",
    "create_point_defect",
    "create_vacancy",
    "create_substitution",
    "create_interstitial",
]


def create_point_defect(
    host_material: Material,
    defect_type: PointDefectTypeEnum,
    coordinate: Optional[List[float]] = None,
    resolution_method: AtomPlacementMethodEnum = AtomPlacementMethodEnum.COORDINATE,
    element: Optional[str] = None,
    builder_parameters: Optional[MergeBuilderParameters] = None,
) -> Material:
    """
    Create a point defect in a material using the new merge-based infrastructure.

    Args:
        host_material: The host crystal material
        defect_type: Type of defect (vacancy, substitution, interstitial)
        coordinate: Optional coordinate for the defect
        resolution_method: Method to resolve coordinates
        element: Element for substitution/interstitial defects
        builder_parameters: Optional builder parameters

    Returns:
        Material with the point defect applied
    """

    from mat3ra.made.tools.analyze.point_defect import PointDefectMaterialAnalyzer

    analyzer = PointDefectMaterialAnalyzer()
    configuration = analyzer.create_point_defect_configuration(
        host_material=host_material,
        defect_type=defect_type,
        coordinate=coordinate,
        resolution_method=resolution_method,
        element=element,
    )

    builder = PointDefectBuilder(build_parameters=builder_parameters)
    return builder.get_material(configuration)


def create_vacancy(
    host_material: Material,
    coordinate: Optional[List[float]] = None,
    resolution_method: AtomPlacementMethodEnum = AtomPlacementMethodEnum.COORDINATE,
    builder_parameters: Optional[MergeBuilderParameters] = None,
) -> Material:
    """
    Create a vacancy defect in a material.

    Args:
        host_material: The host crystal material
        coordinate: Optional coordinate for the vacancy
        resolution_method: Method to resolve coordinates
        builder_parameters: Optional builder parameters

    Returns:
        Material with the vacancy defect
    """

    return create_point_defect(
        host_material=host_material,
        defect_type=PointDefectTypeEnum.VACANCY,
        coordinate=coordinate,
        resolution_method=resolution_method,
        element=None,
        builder_parameters=builder_parameters,
    )


def create_substitution(
    host_material: Material,
    element: str,
    coordinate: Optional[List[float]] = None,
    resolution_method: AtomPlacementMethodEnum = AtomPlacementMethodEnum.COORDINATE,
    builder_parameters: Optional[MergeBuilderParameters] = None,
) -> Material:
    """
    Create a substitution defect in a material.

    Args:
        host_material: The host crystal material
        element: The element to substitute with
        coordinate: Optional coordinate for the substitution
        resolution_method: Method to resolve coordinates
        builder_parameters: Optional builder parameters

    Returns:
        Material with the substitution defect
    """

    return create_point_defect(
        host_material=host_material,
        defect_type=PointDefectTypeEnum.SUBSTITUTION,
        coordinate=coordinate,
        resolution_method=resolution_method,
        element=element,
        builder_parameters=builder_parameters,
    )


def create_interstitial(
    host_material: Material,
    element: str,
    coordinate: Optional[List[float]] = None,
    resolution_method: AtomPlacementMethodEnum = AtomPlacementMethodEnum.COORDINATE,
    builder_parameters: Optional[MergeBuilderParameters] = None,
) -> Material:
    """
    Create an interstitial defect in a material.

    Args:
        host_material: The host crystal material
        element: Element to add as interstitial
        coordinate: Optional coordinate for the interstitial
        resolution_method: Method to resolve coordinates
        builder_parameters: Optional builder parameters

    Returns:
        Material with the interstitial defect
    """

    return create_point_defect(
        host_material=host_material,
        defect_type=PointDefectTypeEnum.INTERSTITIAL,
        coordinate=coordinate,
        resolution_method=resolution_method,
        element=element,
        builder_parameters=builder_parameters,
    )
