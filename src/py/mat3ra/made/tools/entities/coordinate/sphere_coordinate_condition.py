from typing import List

from mat3ra.esse.models.core.reusable.coordinate_conditions import SphereCoordinateConditionSchema
from pydantic import Field

from .coordinate_condition import CoordinateCondition
from .coordinate_functions import is_coordinate_in_sphere


class SphereCoordinateCondition(SphereCoordinateConditionSchema, CoordinateCondition):
    center_coordinate: List[float] = Field(default_factory=lambda: [0.5, 0.5, 0.5])
    radius: float = 0.25

    def condition(self, coordinate: List[float]) -> bool:
        return is_coordinate_in_sphere(coordinate, self.center_coordinate, self.radius)
