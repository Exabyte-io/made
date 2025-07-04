from typing import List

from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.point_defect import PointDefectAnalyzer
from mat3ra.made.tools.build.defect.enums import AtomPlacementMethodEnum, PointDefectTypeEnum
from mat3ra.made.tools.build.defect.point.builders import (
    VacancyDefectBuilder,
    SubstitutionalDefectBuilder,
    InterstitialDefectBuilder,
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
    analyzer = PointDefectAnalyzer(material=material, resolution_method=placement_method)
    config = analyzer.get_configuration(defect_type=PointDefectTypeEnum.VACANCY, coordinate=coordinate)
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
    analyzer = PointDefectAnalyzer(material=material, resolution_method=placement_method)
    config = analyzer.get_configuration(
        defect_type=PointDefectTypeEnum.SUBSTITUTION, coordinate=coordinate, element=element
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
    analyzer = PointDefectAnalyzer(material=material, resolution_method=placement_method)
    config = analyzer.get_configuration(
        defect_type=PointDefectTypeEnum.INTERSTITIAL, coordinate=coordinate, element=element
    )
    builder = InterstitialDefectBuilder()
    return builder.get_material(config)
