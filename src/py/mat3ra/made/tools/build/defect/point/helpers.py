from types import SimpleNamespace
from typing import List, Union

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
from ... import MaterialWithBuildMetadata
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


def create_point_defect_vacancy(
    material: Union[Material, MaterialWithBuildMetadata],
    coordinate: List[float],
    placement_method: str,
    use_cartesian_coordinates: bool = False,
) -> Material:
    """
    Create a vacancy defect in the given material at the specified coordinate.

    Args:
        material (Material): The host material.
        coordinate (List[float]): The coordinate where the vacancy will be created.
        placement_method (VacancyPlacementMethodEnum): Method to resolve the final coordinate.
        use_cartesian_coordinates (bool): Whether the input coordinate is in Cartesian units.

    Returns:
        Material: A new material with the vacancy defect.
    """
    if placement_method not in [e.value for e in VacancyPlacementMethodEnum]:
        raise ValueError(f"Unsupported placement method: {placement_method}")

    # Convert coordinate to crystal if needed
    if use_cartesian_coordinates:
        coordinate = material.basis.cell.convert_point_to_crystal(coordinate)

    analyzer = CrystalSiteAnalyzer(material=material, coordinate=coordinate)
    resolved_coordinate = analyzer.closest_site_coordinate
    config = VacancyDefectConfiguration.from_parameters(crystal=material, coordinate=resolved_coordinate)
    builder = VacancyDefectBuilder()
    return builder.get_material(config)


def create_point_defect_substitution(
    material: Union[Material, MaterialWithBuildMetadata],
    coordinate: List[float],
    element: str,
    placement_method: str,
    use_cartesian_coordinates: bool = False,
) -> Material:
    """
    Create a substitution defect in the given material.

    Args:
        material (Material): The host material.
        coordinate (List[float]): The coordinate of the atom to be substituted.
        element (str): The chemical element to substitute with.
        placement_method (SubstitutionPlacementMethodEnum): Method to resolve the final coordinate.
        use_cartesian_coordinates (bool): Whether the input coordinate is in Cartesian units.

    Returns:
        Material: A new material with the substitution defect.
    """
    if placement_method not in [e.value for e in SubstitutionPlacementMethodEnum]:
        raise ValueError(f"Unsupported placement method: {placement_method}")

    if use_cartesian_coordinates:
        coordinate = material.basis.cell.convert_point_to_crystal(coordinate)

    analyzer = CrystalSiteAnalyzer(material=material, coordinate=coordinate)
    resolved_coordinate = analyzer.closest_site_coordinate
    config = SubstitutionalDefectConfiguration.from_parameters(
        crystal=material, coordinate=resolved_coordinate, element=element
    )
    builder = SubstitutionalDefectBuilder()
    return builder.get_material(config)


def create_point_defect_interstitial(
    material: Union[Material, MaterialWithBuildMetadata],
    coordinate: List[float],
    element: str,
    placement_method: str,
    use_cartesian_coordinates: bool = False,
) -> Material:
    """
    Create an interstitial defect in the given material.

    Args:
        material (Material): The host material.
        coordinate (List[float]): The coordinate where the interstitial atom will be placed.
        element (str): The chemical element of the interstitial atom.
        placement_method (InterstitialPlacementMethodEnum): Method to resolve the final coordinate.
        use_cartesian_coordinates (bool): Whether the input coordinate is in Cartesian units.

    Returns:
        Material: A new material with the interstitial defect.
    """
    if use_cartesian_coordinates:
        coordinate = material.basis.cell.convert_point_to_crystal(coordinate)

    if placement_method == InterstitialPlacementMethodEnum.VORONOI_SITE.value:
        analyzer = VoronoiCrystalSiteAnalyzer(material=material, coordinate=coordinate)
        resolved_coordinate = analyzer.voronoi_site_coordinate
    elif placement_method == InterstitialPlacementMethodEnum.EXACT_COORDINATE.value:
        analyzer = CrystalSiteAnalyzer(material=material, coordinate=coordinate)
        resolved_coordinate = analyzer.exact_coordinate
    else:
        raise ValueError(f"Unsupported placement method for interstitial: {placement_method}")

    config = InterstitialDefectConfiguration.from_parameters(
        crystal=material, coordinate=resolved_coordinate, element=element
    )
    builder = InterstitialDefectBuilder()
    return builder.get_material(config)


def create_multiple_defects(
    material: Union[Material, MaterialWithBuildMetadata],
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
            - use_cartesian_coordinates: bool (optional, defaults to False)
                Whether the input coordinates are in Cartesian units.

    Returns:
        Material: A new material with all defects applied.
    """
    current_material = material

    for defect_dict in defect_dicts:
        defect_configuration = SimpleNamespace(**defect_dict)
        defect_type = defect_configuration.type

        if defect_type not in [e.value for e in PointDefectTypeEnum]:
            raise ValueError(f"Unsupported defect type: {defect_configuration.type}")

        use_cartesian = getattr(defect_configuration, "use_cartesian_coordinates", False)

        if defect_type == "vacancy":
            current_material = create_point_defect_vacancy(
                current_material,
                coordinate=defect_configuration.coordinate,
                placement_method=defect_configuration.placement_method or VacancyPlacementMethodEnum.CLOSEST_SITE.value,
                use_cartesian_coordinates=use_cartesian,
            )

        elif defect_type == "substitution":
            current_material = create_point_defect_substitution(
                current_material,
                coordinate=defect_configuration.coordinate,
                element=defect_configuration.element,
                placement_method=defect_configuration.placement_method
                or SubstitutionPlacementMethodEnum.CLOSEST_SITE.value,
                use_cartesian_coordinates=use_cartesian,
            )

        elif defect_type == "interstitial":
            current_material = create_point_defect_interstitial(
                current_material,
                coordinate=defect_configuration.coordinate,
                element=defect_configuration.element,
                placement_method=defect_configuration.placement_method
                or InterstitialPlacementMethodEnum.EXACT_COORDINATE.value,
                use_cartesian_coordinates=use_cartesian,
            )

    return current_material
