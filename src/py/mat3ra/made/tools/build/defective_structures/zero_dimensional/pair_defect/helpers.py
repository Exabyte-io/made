from typing import List, Optional, Union

from mat3ra.made.material import Material
from .builder import PairDefectBuilder
from .configuration import PairDefectConfiguration
from ..point_defect import (
    PointDefectTypeEnum,
    VacancyPlacementMethodEnum,
    SubstitutionPlacementMethodEnum,
    InterstitialPlacementMethodEnum,
    AtomPlacementMethodEnum,
    create_defect_configuration,
)
from .....build_components import MaterialWithBuildMetadata


def create_defect_pair(
    material: Union[Material, MaterialWithBuildMetadata],
    defect_type_1: Optional[PointDefectTypeEnum] = None,
    coordinate_1: Optional[List[float]] = None,
    element_1: Optional[str] = None,
    placement_method_1: Union[
        VacancyPlacementMethodEnum,
        SubstitutionPlacementMethodEnum,
        InterstitialPlacementMethodEnum,
        AtomPlacementMethodEnum,
    ] = AtomPlacementMethodEnum.EXACT_COORDINATE,
    defect_type_2: Optional[PointDefectTypeEnum] = None,
    coordinate_2: Optional[List[float]] = None,
    element_2: Optional[str] = None,
    placement_method_2: Union[
        VacancyPlacementMethodEnum,
        SubstitutionPlacementMethodEnum,
        InterstitialPlacementMethodEnum,
        AtomPlacementMethodEnum,
    ] = AtomPlacementMethodEnum.EXACT_COORDINATE,
    use_cartesian_coordinates: bool = False,
) -> Material:
    """
    Create a pair defect in the given material.

    Args:
        material: The host material.
        defect_type_1: Type of the first defect.
        coordinate_1: Coordinate for the first defect.
        element_1: Element for substitution/interstitial defects.
        placement_method_1: Method to resolve the final coordinate for the first defect.
        defect_type_2: Type of the second defect.
        coordinate_2: Coordinate for the second defect.
        element_2: Element for substitution/interstitial defects.
        placement_method_2: Method to resolve the final coordinate for the second defect.
        use_cartesian_coordinates: Whether the coordinates are in Cartesian units.

    Returns:
        material: Union[Material, MaterialWithBuildMetadata] with the pair defect applied.
    """
    if not defect_type_1 or not coordinate_1:
        raise ValueError("defect_type_1 and coordinate_1 must be provided.")
    configuration_1 = create_defect_configuration(
        material,
        defect_type_1,
        coordinate_1,
        element_1,
        placement_method_1,
        use_cartesian_coordinates,
    )
    if not defect_type_2 or not coordinate_2:
        raise ValueError("defect_type_2 and coordinate_2 must be provided.")
    configuration_2 = create_defect_configuration(
        material, defect_type_2, coordinate_2, element_2, placement_method_2, use_cartesian_coordinates
    )

    pair_config = PairDefectConfiguration.from_parameters(
        crystal=material,
        primary_defect_configuration=configuration_1,
        secondary_defect_configuration=configuration_2,
    )

    builder = PairDefectBuilder()
    return builder.get_material(pair_config)
