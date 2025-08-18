from typing import Union, List

from mat3ra.made.material import Material
from .builder import SubstitutionalDefectBuilder
from .configuration import SubstitutionalDefectConfiguration
from .substitution_placement_method_enum import SubstitutionPlacementMethodEnum
from ......analyze.crystal_site.crystal_site_analyzer import CrystalSiteAnalyzer
from ......build_components import MaterialWithBuildMetadata


def create_defect_point_substitution(
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
