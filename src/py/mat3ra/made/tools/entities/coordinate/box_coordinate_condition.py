from typing import List

from mat3ra.esse.models.core.reusable.coordinate_conditions import BoxCoordinateConditionSchema
from pydantic import Field

from .coordinate_condition import CoordinateCondition
from .coordinate_functions import is_coordinate_in_box


class BoxCoordinateCondition(BoxCoordinateConditionSchema, CoordinateCondition):
    min_coordinate: List[float] = Field(default_factory=lambda: [0, 0, 0])
    max_coordinate: List[float] = Field(default_factory=lambda: [1, 1, 1])

    def condition(self, coordinate: List[float]) -> bool:
        return is_coordinate_in_box(coordinate, self.min_coordinate, self.max_coordinate)
