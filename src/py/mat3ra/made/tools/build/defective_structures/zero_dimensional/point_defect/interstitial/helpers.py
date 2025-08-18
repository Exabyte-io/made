from typing import Union, List

from mat3ra.made.material import Material
from ......analyze.crystal_site.crystal_site_analyzer import CrystalSiteAnalyzer
from ......analyze.crystal_site.voronoi_crystal_site_analyzer import VoronoiCrystalSiteAnalyzer
from .builder import InterstitialDefectBuilder
from .configuration import InterstitialDefectConfiguration
from .interstitial_placement_method_enum import InterstitialPlacementMethodEnum
from ......build_components import MaterialWithBuildMetadata


def create_defect_point_interstitial(
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
