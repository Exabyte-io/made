from typing import List

from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.point_defect import CrystalSiteAnalyzer, VoronoiCrystalSiteAnalyzer
from mat3ra.made.tools.build.defect.enums import AtomPlacementMethodEnum, PointDefectTypeEnum
from mat3ra.made.tools.build.defect.point.builders import (
    VacancyDefectBuilder,
    SubstitutionalDefectBuilder,
    InterstitialDefectBuilder,
)
from mat3ra.made.tools.build.defect.point.configuration import (
    VacancyDefectConfiguration,
    SubstitutionalDefectConfiguration,
    InterstitialDefectConfiguration,
)


def create_vacancy_defect(
    material: Material,
    coordinate: List[float],
    placement_method: AtomPlacementMethodEnum = AtomPlacementMethodEnum.COORDINATE,
) -> Material:
    """
    Create a vacancy defect in the given material at the specified coordinate.

    Args:
        material (Material): The host material.
        coordinate (List[float]): The coordinate where the vacancy will be created.
        placement_method (AtomPlacementMethodEnum): Method to resolve the final coordinate.

    Returns:
        Material: A new material with the vacancy defect.
    """
    if placement_method == AtomPlacementMethodEnum.VORONOI_SITE:
        analyzer = VoronoiCrystalSiteAnalyzer(material=material, coordinate=coordinate)
    else:
        analyzer = CrystalSiteAnalyzer(material=material, coordinate=coordinate)
    
    resolved_coordinate = _get_resolved_crystal_site_coordinate(analyzer, placement_method)
    config = VacancyDefectConfiguration.from_parameters(crystal=material, coordinate=resolved_coordinate)
    builder = VacancyDefectBuilder()
    return builder.get_material(config)


def create_substitution_defect(
    material: Material,
    coordinate: List[float],
    element: str,
    placement_method: AtomPlacementMethodEnum = AtomPlacementMethodEnum.COORDINATE,
) -> Material:
    """
    Create a substitution defect in the given material.

    Args:
        material (Material): The host material.
        coordinate (List[float]): The coordinate of the atom to be substituted.
        element (str): The chemical element to substitute with.
        placement_method (AtomPlacementMethodEnum): Method to resolve the final coordinate.

    Returns:
        Material: A new material with the substitution defect.
    """
    if placement_method == AtomPlacementMethodEnum.VORONOI_SITE:
        analyzer = VoronoiCrystalSiteAnalyzer(material=material, coordinate=coordinate)
    else:
        analyzer = CrystalSiteAnalyzer(material=material, coordinate=coordinate)
    
    resolved_coordinate = _get_resolved_crystal_site_coordinate(analyzer, placement_method)
    config = SubstitutionalDefectConfiguration.from_parameters(
        crystal=material, coordinate=resolved_coordinate, element=element
    )
    builder = SubstitutionalDefectBuilder()
    return builder.get_material(config)


def create_interstitial_defect(
    material: Material,
    coordinate: List[float],
    element: str,
    placement_method: AtomPlacementMethodEnum = AtomPlacementMethodEnum.COORDINATE,
) -> Material:
    """
    Create an interstitial defect in the given material.

    Args:
        material (Material): The host material.
        coordinate (List[float]): The coordinate where the interstitial atom will be placed.
        element (str): The chemical element of the interstitial atom.
        placement_method (AtomPlacementMethodEnum): Method to resolve the final coordinate.

    Returns:
        Material: A new material with the interstitial defect.
    """
    if placement_method == AtomPlacementMethodEnum.VORONOI_SITE:
        analyzer = VoronoiCrystalSiteAnalyzer(material=material, coordinate=coordinate)
    else:
        analyzer = CrystalSiteAnalyzer(material=material, coordinate=coordinate)
    
    resolved_coordinate = _get_resolved_crystal_site_coordinate(analyzer, placement_method)
    config = InterstitialDefectConfiguration.from_parameters(
        crystal=material, coordinate=resolved_coordinate, element=element
    )
    builder = InterstitialDefectBuilder()
    return builder.get_material(config)


def _get_resolved_crystal_site_coordinate(analyzer, placement_method: AtomPlacementMethodEnum) -> List[float]:
    """Get the resolved coordinate based on the placement method."""
    if placement_method == AtomPlacementMethodEnum.COORDINATE:
        return analyzer.coordinate_resolution
    elif placement_method == AtomPlacementMethodEnum.CLOSEST_SITE:
        return analyzer.closest_site_resolution
    elif placement_method == AtomPlacementMethodEnum.NEW_CRYSTAL_SITE:
        return analyzer.new_crystal_site_resolution
    elif placement_method == AtomPlacementMethodEnum.EQUIDISTANT:
        return analyzer.equidistant_resolution
    elif placement_method == AtomPlacementMethodEnum.VORONOI_SITE:
        return analyzer.voronoi_site_resolution
    else:
        raise ValueError(f"Unknown atom placement method: {placement_method}")
