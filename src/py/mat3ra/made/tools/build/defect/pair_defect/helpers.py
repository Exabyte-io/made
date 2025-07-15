from typing import List

from mat3ra.made.material import Material
from .builders import PairDefectBuilder
from .configuration import PairDefectConfiguration
from ..enums import PointDefectTypeEnum
from ..factories import create_defect_configuration


def create_pair_defect(
    material: Material,
    defect_type_1: PointDefectTypeEnum = None,
    coordinate_1: List[float] = None,
    element_1: str = None,
    defect_type_2: PointDefectTypeEnum = None,
    coordinate_2: List[float] = None,
    element_2: str = None,
) -> Material:
    """
    Create a pair defect in the given material.

    Args:
        material: The host material.
        defect_type_1: Type of the first defect.
        coordinate_1: Coordinate for the first defect.
        element_1: Element for substitution/interstitial defects.
        defect_type_2: Type of the second defect.
        coordinate_2: Coordinate for the second defect.
        element_2: Element for substitution/interstitial defects.

    Returns:
        Material: Material with the pair defect applied.
    """
    configuration_1 = create_defect_configuration(material, defect_type_1, coordinate_1, element_1)
    configuration_2 = create_defect_configuration(material, defect_type_2, coordinate_2, element_2)

    pair_config = PairDefectConfiguration.from_parameters(
        crystal=material,
        primary_defect_configuration=configuration_1,
        secondary_defect_configuration=configuration_2,
    )

    builder = PairDefectBuilder()
    return builder.get_material(pair_config)
