from typing import List

from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.crystal_site import CrystalSiteAnalyzer
from mat3ra.made.tools.build.defect.point.builders import PointDefectBuilder
from mat3ra.made.tools.build.defect.point.configuration import PointDefectConfiguration


def create_multiple_defects(
    material: Material,
    defect_configurations: List[PointDefectConfiguration],
    coord_resolution_method: str = "closest_site",
) -> Material:
    """
    Create a material with multiple defects.

    Args:
        material: The host material.
        defect_configurations: List of point defect configurations.
        coord_resolution_method: Method to resolve coordinates ("closest_site", "exact_coordinate", etc.).

    Returns:
        Material: Material with multiple defects applied.
    """
    current_material = material

    for defect_config in defect_configurations:
        builder = PointDefectBuilder()
        current_material = builder.get_material(defect_config)

    return current_material


def resolve_coordinate_by_method(
    material: Material,
    coordinate: List[float],
    method: str,
) -> List[float]:
    """
    Resolve coordinate using the specified method.

    Args:
        material: The host material.
        coordinate: The input coordinate.
        method: The resolution method.

    Returns:
        List[float]: The resolved coordinate.
    """
    analyzer = CrystalSiteAnalyzer(material=material, coordinate=coordinate)

    if method == "closest_site":
        return analyzer.closest_site_coordinate
    elif method == "exact_coordinate":
        return analyzer.exact_coordinate
    elif method == "equidistant":
        return analyzer.get_equidistant_coordinate()
    else:
        raise ValueError(f"Unknown coordinate resolution method: {method}")
