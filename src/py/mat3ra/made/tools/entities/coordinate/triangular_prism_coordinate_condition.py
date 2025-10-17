from typing import List

from mat3ra.esse.models.core.reusable.coordinate_conditions.triangular_prism import (
    TriangularPrismCoordinateConditionSchema,
)

from .coordinate_condition import CoordinateCondition
from .coordinate_functions import is_coordinate_in_triangular_prism


class TriangularPrismCoordinateCondition(TriangularPrismCoordinateConditionSchema, CoordinateCondition):
    position_on_surface_1: List[float] = [0, 0]
    position_on_surface_2: List[float] = [1, 0]
    position_on_surface_3: List[float] = [0, 1]
    min_z: float = 0
    max_z: float = 1

    def condition(self, coordinate: List[float]) -> bool:
        return is_coordinate_in_triangular_prism(
            coordinate,
            self.position_on_surface_1,
            self.position_on_surface_2,
            self.position_on_surface_3,
            self.min_z,
            self.max_z,
        )
