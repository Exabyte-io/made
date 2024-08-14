from typing import List

from pydantic import BaseModel

AXIS_TO_INDEX_MAP = {"x": 0, "y": 1, "z": 2}
EQUATION_RANGE_COEFFICIENT = 5


class FunctionHolder(BaseModel):
    def apply_function(self, coordinate: List[float]) -> float:
        """
        Get the value of the function at the given coordinate.
        """
        raise NotImplementedError

    def calculate_derivative(self, coordinate: List[float]) -> float:
        """
        Get the derivative of the function at the given coordinate
        """
        raise NotImplementedError

    def calculate_arc_length(self, a: float, b: float) -> float:
        """
        Get the arc length of the function between a and b.
        """
        raise NotImplementedError

    def get_json(self) -> dict:
        """
        Get the json representation of the function holder.
        """
        raise NotImplementedError
