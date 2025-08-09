from typing import List

from mat3ra.esse.models.core.reusable.coordinate_conditions import CylinderCoordinateConditionSchema
from pydantic import Field

from .coordinate_condition import CoordinateCondition
from .coordinate_functions import is_coordinate_in_cylinder


class CylinderCoordinateCondition(CylinderCoordinateConditionSchema, CoordinateCondition):
    center_position: List[float] = Field(default_factory=lambda: [0.5, 0.5])
    radius: float = 0.25
    min_z: float = 0
    max_z: float = 1

    def condition(self, coordinate: List[float]) -> bool:
        return is_coordinate_in_cylinder(coordinate, self.center_position, self.radius, self.min_z, self.max_z)
