from typing import List, Literal, TypedDict
from typing_extensions import NotRequired


class PointDefectDict(TypedDict):
    """
    TypedDict for point defect configurations used with create_multiple_defects.

    Required fields:
        type: The type of defect ("vacancy", "substitution", "interstitial")
        coordinate: Position coordinates as [x, y, z] list

    Optional fields:
        element: Chemical element (required for substitution and interstitial)
        placement_method: Method for placing the defect
        use_cartesian_coordinates: Whether coordinates are in Cartesian units
    """

    type: Literal["vacancy", "substitution", "interstitial"]
    coordinate: List[float]
    element: NotRequired[str]  # Required for substitution and interstitial
    placement_method: NotRequired[str]  # "CLOSEST_SITE", "EXACT_COORDINATE", "VORONOI_SITE"
    use_cartesian_coordinates: NotRequired[bool]  # defaults to False
