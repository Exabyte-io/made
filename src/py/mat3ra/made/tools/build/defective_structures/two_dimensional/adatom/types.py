from typing import List, TypedDict
from typing_extensions import NotRequired


class AdatomDefectDict(TypedDict):
    """
    TypedDict for adatom defect configurations used with create_multiple_adatom_defects.

    Required fields:
        element: Chemical element for the adatom
        coordinate: Position on surface as [x, y] list
        distance_z: Distance above surface in Angstroms

    Optional fields:
        use_cartesian_coordinates: Whether coordinates are in Cartesian units
    """

    element: str
    coordinate: List[float]  # [x, y] position on surface
    distance_z: float
    use_cartesian_coordinates: NotRequired[bool]  # defaults to False
