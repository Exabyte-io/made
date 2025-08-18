from typing import List, Union

from mat3ra.made.material import Material
from .interstitial.helpers import create_defect_point_interstitial
from .interstitial.interstitial_placement_method_enum import InterstitialPlacementMethodEnum
from .point_defect_type_enum import PointDefectTypeEnum
from .substitutional.helpers import create_defect_point_substitution
from .substitutional.substitution_placement_method_enum import SubstitutionPlacementMethodEnum
from .types import PointDefectDict
from .vacancy.helpers import create_defect_point_vacancy
from .vacancy.vacancy_placement_method_enum import VacancyPlacementMethodEnum
from .....build_components import MaterialWithBuildMetadata

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


def create_multiple_defects(
    material: Union[Material, MaterialWithBuildMetadata],
    defect_dicts: List[PointDefectDict],
) -> Material:
    """
    Create multiple point defects from a list of PointDefectDict.

    Args:
        material (Material): The host material.
        defect_dicts (List[PointDefectDict]): List of defect dictionaries with keys:
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
        defect_type = defect_dict.type

        if defect_type not in [e.value for e in PointDefectTypeEnum]:
            raise ValueError(f"Unsupported defect type: {defect_dict.type}")

        use_cartesian = getattr(defect_dict, "use_cartesian_coordinates", False)

        if defect_type == "vacancy":
            current_material = create_defect_point_vacancy(
                current_material,
                coordinate=defect_dict.coordinate,
                placement_method=defect_dict.placement_method or VacancyPlacementMethodEnum.CLOSEST_SITE.value,
                use_cartesian_coordinates=use_cartesian,
            )

        elif defect_type == "substitution":
            current_material = create_defect_point_substitution(
                current_material,
                coordinate=defect_dict.coordinate,
                element=defect_dict.element,  # type: ignore # Pydantic validates this at creation
                placement_method=defect_dict.placement_method or SubstitutionPlacementMethodEnum.CLOSEST_SITE.value,
                use_cartesian_coordinates=use_cartesian,
            )

        elif defect_type == "interstitial":
            current_material = create_defect_point_interstitial(
                current_material,
                coordinate=defect_dict.coordinate,
                element=defect_dict.element,  # type: ignore # Pydantic validates this at creation
                placement_method=defect_dict.placement_method or InterstitialPlacementMethodEnum.EXACT_COORDINATE.value,
                use_cartesian_coordinates=use_cartesian,
            )

    return current_material
