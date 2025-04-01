from typing import Literal

import numpy as np
from mat3ra.code.array_with_ids import RoundedArrayWithIds
from mat3ra.code.value_with_id import RoundedValueWithId
from mat3ra.esse.models.core.abstract.point import PointSchema


class Point(PointSchema):
    pass


class Coordinate(RoundedValueWithId):
    value: Point

    def get_value_along_axis(self, axis: Literal["x", "y", "z"] = "z"):
        return self.value.root[{"x": 0, "y": 1, "z": 2}[axis]]


class Coordinates(RoundedArrayWithIds):
    def get_values_along_axis(
        self,
        axis: Literal["x", "y", "z"] = "z",
    ):
        values_along_axis = [Coordinate(value=coord).get_value_along_axis(axis) for coord in self.values]
        return values_along_axis

    def get_max_value_along_axis(
        self,
        axis: Literal["x", "y", "z"] = "z",
    ):
        return np.max(self.get_values_along_axis(axis))

    def get_min_value_along_axis(
        self,
        axis: Literal["x", "y", "z"] = "z",
    ):
        return np.min(self.get_values_along_axis(axis))

    def get_extremum_value_along_axis(
        self,
        extremum: Literal["max", "min"] = "max",
        axis: Literal["x", "y", "z"] = "z",
    ):
        if extremum == "max":
            return self.get_max_value_along_axis(axis)
        return self.get_min_value_along_axis(axis)
