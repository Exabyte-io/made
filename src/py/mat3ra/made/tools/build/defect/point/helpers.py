from typing import List

from mat3ra.made.material import Material
from .builders import (
    VacancyDefectBuilder,
    SubstitutionalDefectBuilder,
    InterstitialDefectBuilder,
)
from .configuration import (
    VacancyDefectConfiguration,
    SubstitutionalDefectConfiguration,
    InterstitialDefectConfiguration,
)
from ..enums import (
    PointDefectTypeEnum,
    VacancyPlacementMethodEnum,
    SubstitutionPlacementMethodEnum,
    InterstitialPlacementMethodEnum,
)
from ....analyze.crystal_site import CrystalSiteAnalyzer, VoronoiCrystalSiteAnalyzer


def create_point_defect_vacancy(
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


def create_point_defect_substitution(
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


def create_point_defect_interstitial(
    material: Material,
    coordinate: List[float],
    element: str,
    placement_method: InterstitialPlacementMethodEnum = InterstitialPlacementMethodEnum.EXACT_COORDINATE,
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
    if placement_method == InterstitialPlacementMethodEnum.VORONOI_SITE:
        analyzer = VoronoiCrystalSiteAnalyzer(material=material, coordinate=coordinate)
        resolved_coordinate = analyzer.voronoi_site_coordinate
    else:
        analyzer = CrystalSiteAnalyzer(material=material, coordinate=coordinate)
        resolved_coordinate = analyzer.exact_coordinate

    config = InterstitialDefectConfiguration.from_parameters(
        crystal=material, coordinate=resolved_coordinate, element=element
    )
    builder = InterstitialDefectBuilder()
    return builder.get_material(config)


def create_multiple_defects(
    material: Material,
    defect_dicts: List[dict],
) -> Material:
    """
    Create multiple point defects from a list of dictionaries.

    Args:
        material (Material): The host material.
        defect_dicts (List[dict]): List of defect dictionaries with keys:
            - type: str ("vacancy", "substitution", "interstitial")
            - coordinate: List[float]
            - element: str (required for substitution and interstitial)
            - resolution_method: str (placement method)

    Returns:
        Material: A new material with all defects applied.
    """
    current_material = material

    for defect_dict in defect_dicts:
        defect_type = defect_dict["type"]
        coordinate = defect_dict["coordinate"]
        resolution_method = defect_dict.get("resolution_method")

        if defect_type == PointDefectTypeEnum.VACANCY:
            placement_method = VacancyPlacementMethodEnum.CLOSEST_SITE
            if resolution_method:
                placement_method = VacancyPlacementMethodEnum(resolution_method)
            current_material = create_point_defect_vacancy(current_material, coordinate, placement_method)

        elif defect_type == PointDefectTypeEnum.SUBSTITUTION:
            element = defect_dict["element"]
            placement_method = SubstitutionPlacementMethodEnum.CLOSEST_SITE
            if resolution_method:
                placement_method = SubstitutionPlacementMethodEnum(resolution_method)
            current_material = create_point_defect_substitution(current_material, coordinate, element, placement_method)

        elif defect_type == PointDefectTypeEnum.INTERSTITIAL:
            element = defect_dict["element"]
            placement_method = InterstitialPlacementMethodEnum.EXACT_COORDINATE
            if resolution_method:
                placement_method = InterstitialPlacementMethodEnum(resolution_method)
            current_material = create_point_defect_interstitial(current_material, coordinate, element, placement_method)

        else:
            raise ValueError(f"Unsupported defect type: {defect_type}")

    return current_material
