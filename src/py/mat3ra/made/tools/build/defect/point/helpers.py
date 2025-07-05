from typing import List

from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.crystal_site import CrystalSiteAnalyzer, VoronoiCrystalSiteAnalyzer
from mat3ra.made.tools.build.defect.enums import (
    VacancyPlacementMethodEnum,
    SubstitutionPlacementMethodEnum,
    InterstitialPlacementMethodEnum,
)
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
    placement_method: VacancyPlacementMethodEnum = VacancyPlacementMethodEnum.CLOSEST_SITE,
) -> Material:
    """
    Create a vacancy defect in the given material at the specified coordinate.

    Args:
        material (Material): The host material.
        coordinate (List[float]): The coordinate where the vacancy will be created.
        placement_method (VacancyPlacementMethodEnum): Method to resolve the final coordinate.

    Returns:
        Material: A new material with the vacancy defect.
    """
    analyzer = CrystalSiteAnalyzer(material=material, coordinate=coordinate)
    if placement_method == VacancyPlacementMethodEnum.CLOSEST_SITE:
        resolved_coordinate = analyzer.closest_site_coordinate
    else:
        raise NotImplementedError(f"Vacancy placement method '{placement_method}' is not implemented.")
    config = VacancyDefectConfiguration.from_parameters(crystal=material, coordinate=resolved_coordinate)
    builder = VacancyDefectBuilder()
    return builder.get_material(config)


def create_substitution_defect(
    material: Material,
    coordinate: List[float],
    element: str,
    placement_method: SubstitutionPlacementMethodEnum = SubstitutionPlacementMethodEnum.CLOSEST_SITE,
) -> Material:
    """
    Create a substitution defect in the given material.

    Args:
        material (Material): The host material.
        coordinate (List[float]): The coordinate of the atom to be substituted.
        element (str): The chemical element to substitute with.
        placement_method (SubstitutionPlacementMethodEnum): Method to resolve the final coordinate.

    Returns:
        Material: A new material with the substitution defect.
    """
    analyzer = CrystalSiteAnalyzer(material=material, coordinate=coordinate)
    if placement_method == SubstitutionPlacementMethodEnum.CLOSEST_SITE:
        resolved_coordinate = analyzer.closest_site_coordinate
    else:
        raise NotImplementedError(f"Substitution placement method '{placement_method}' is not implemented.")
    config = SubstitutionalDefectConfiguration.from_parameters(
        crystal=material, coordinate=resolved_coordinate, element=element
    )
    builder = SubstitutionalDefectBuilder()
    return builder.get_material(config)


def create_interstitial_defect(
    material: Material,
    coordinate: List[float],
    element: str,
    placement_method: InterstitialPlacementMethodEnum = InterstitialPlacementMethodEnum.COORDINATE,
) -> Material:
    """
    Create an interstitial defect in the given material.

    Args:
        material (Material): The host material.
        coordinate (List[float]): The coordinate where the interstitial atom will be placed.
        element (str): The chemical element of the interstitial atom.
        placement_method (InterstitialPlacementMethodEnum): Method to resolve the final coordinate.

    Returns:
        Material: A new material with the interstitial defect.
    """
    if placement_method == InterstitialPlacementMethodEnum.COORDINATE:
        resolved_coordinate = coordinate
    elif placement_method == InterstitialPlacementMethodEnum.VORONOI_SITE:
        analyzer = VoronoiCrystalSiteAnalyzer(material=material, coordinate=coordinate)
        resolved_coordinate = analyzer.voronoi_site_coordinate
    else:
        raise NotImplementedError(f"Interstitial placement method '{placement_method}' is not implemented.")
    config = InterstitialDefectConfiguration.from_parameters(
        crystal=material, coordinate=resolved_coordinate, element=element
    )
    builder = InterstitialDefectBuilder()
    return builder.get_material(config)
