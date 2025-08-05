from typing import List

from .coordinate_condition import CoordinateCondition
from .coordinate_functions import is_coordinate_behind_plane


class PlaneCoordinateCondition(CoordinateCondition):
    plane_normal: List[float]
    plane_point_coordinate: List[float]

    def condition(self, coordinate: List[float]) -> bool:
        return is_coordinate_behind_plane(coordinate, self.plane_normal, self.plane_point_coordinate)
