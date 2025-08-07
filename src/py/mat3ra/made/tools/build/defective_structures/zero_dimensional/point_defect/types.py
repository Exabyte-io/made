from typing import List, Literal, Optional

from pydantic import BaseModel, Field


class PointDefectDict(BaseModel):
    """
    Pydantic model for point defect configurations used with create_multiple_defects.

    Required fields:
        type: The type of defect ("vacancy", "substitution", "interstitial")
        coordinate: Position coordinates as [x, y, z] list

    Optional fields:
        element: Chemical element (required for substitution and interstitial)
        placement_method: Method for placing the defect
        use_cartesian_coordinates: Whether coordinates are in Cartesian units
    """

    type: Literal["vacancy", "substitution", "interstitial"]
    coordinate: List[float] = Field(..., min_items=3, max_items=3, description="Position coordinates as [x, y, z]")
    element: Optional[str] = None
    placement_method: Optional[Literal["closest_site", "exact_coordinate", "voronoi_site"]] = None
    use_cartesian_coordinates: bool = False
