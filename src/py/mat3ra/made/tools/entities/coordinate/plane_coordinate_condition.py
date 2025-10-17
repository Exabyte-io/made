from typing import List

from mat3ra.esse.models.core.reusable.coordinate_conditions.plane import PlaneCoordinateConditionSchema

from .coordinate_condition import CoordinateCondition
from .coordinate_functions import is_coordinate_behind_plane


class PlaneCoordinateCondition(PlaneCoordinateConditionSchema, CoordinateCondition):
    plane_normal: List[float]
    plane_point_coordinate: List[float]

    def condition(self, coordinate: List[float]) -> bool:
        return is_coordinate_behind_plane(coordinate, self.plane_normal, self.plane_point_coordinate)
