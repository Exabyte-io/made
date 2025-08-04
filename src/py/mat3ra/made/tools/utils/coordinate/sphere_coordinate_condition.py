from typing import List

from pydantic import Field

from .coordinate_functions import is_coordinate_in_sphere
from .coordinate_condition import CoordinateCondition


class SphereCoordinateCondition(CoordinateCondition):
    center_coordinate: List[float] = Field(default_factory=lambda: [0.5, 0.5, 0.5])
    radius: float = 0.25

    def condition(self, coordinate: List[float]) -> bool:
        return is_coordinate_in_sphere(coordinate, self.center_coordinate, self.radius)
