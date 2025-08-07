from typing import List

from pydantic import BaseModel, Field


class AdatomDefectDict(BaseModel):
    """
    Pydantic model for adatom defect configurations used with create_multiple_adatom_defects.

    Required fields:
        element: Chemical element for the adatom
        coordinate: Position on surface as [x, y] list
        distance_z: Distance above surface in Angstroms

    Optional fields:
        use_cartesian_coordinates: Whether coordinates are in Cartesian units
    """

    element: str = Field(..., description="Chemical element for the adatom")
    coordinate_2d: List[float] = Field(..., min_items=2, max_items=2, description="Position on surface as [x, y]")
    distance_z: float = Field(..., gt=0, description="Distance above surface in Angstroms")
    use_cartesian_coordinates: bool = False
