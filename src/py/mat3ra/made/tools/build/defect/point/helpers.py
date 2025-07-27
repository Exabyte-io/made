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

DEFECT_TYPE_MAPPING = {
    PointDefectTypeEnum.VACANCY: {
        "create_func": "create_point_defect_vacancy",
        "placement_enum": VacancyPlacementMethodEnum,
        "default_method": VacancyPlacementMethodEnum.CLOSEST_SITE,
    },
    PointDefectTypeEnum.SUBSTITUTION: {
        "create_func": "create_point_defect_substitution",
        "placement_enum": SubstitutionPlacementMethodEnum,
        "default_method": SubstitutionPlacementMethodEnum.CLOSEST_SITE,
    },
    PointDefectTypeEnum.INTERSTITIAL: {
        "create_func": "create_point_defect_interstitial",
        "placement_enum": InterstitialPlacementMethodEnum,
        "default_method": InterstitialPlacementMethodEnum.EXACT_COORDINATE,
    },
}


def _convert_placement_method_to_enum(placement_method, placement_enum_class):
    """
    Convert placement method string to enum value.

    Args:
        placement_method: String or enum value
        placement_enum_class: The enum class to convert to

    Returns:
        Enum value
    """
    if isinstance(placement_method, placement_enum_class):
        return placement_method

    if isinstance(placement_method, str):
        try:
            return placement_enum_class[placement_method.upper()]
        except KeyError:
            raise ValueError(f"Invalid placement method '{placement_method}' for {placement_enum_class.__name__}")

    if hasattr(placement_method, "value"):
        for enum_member in placement_enum_class:
            if getattr(enum_member.value, "value", enum_member.value) == placement_method.value:
                return enum_member

    raise ValueError(f"Cannot convert placement method '{placement_method}' to {placement_enum_class.__name__}")


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
            - placement_method: str or enum (optional)
                - For vacancy/substitution: "CLOSEST_SITE"
                - For interstitial:  "EXACT_COORDINATE", "VORONOI_SITE"
                Defaults to "closest_site" for vacancy/substitution and "exact_coordinate" for interstitial.


    Returns:
        Material: A new material with all defects applied.
    """
    current_material = material

    for defect_dict in defect_dicts:
        defect_type = defect_dict["type"]

        if defect_type not in DEFECT_TYPE_MAPPING:
            raise ValueError(f"Unsupported defect type: {defect_type}")

        defect_info = DEFECT_TYPE_MAPPING[defect_type]
        create_func = globals()[defect_info["create_func"]]

        placement_method = defect_dict.get("placement_method") or defect_info["default_method"]

        placement_method_enum = _convert_placement_method_to_enum(placement_method, defect_info["placement_enum"])

        args = [current_material, defect_dict["coordinate"], placement_method_enum]
        if defect_type in (PointDefectTypeEnum.SUBSTITUTION, PointDefectTypeEnum.INTERSTITIAL):
            args.insert(-1, defect_dict["element"])

        current_material = create_func(*args)

    return current_material
