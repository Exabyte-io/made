from typing import List

from mat3ra.code.entity import InMemoryEntityPydantic


class FunctionHolder(InMemoryEntityPydantic):
    def apply_function(self, coordinate: List[float]) -> float:
        """
        Get the value of the function at the given coordinate.
        """
        raise NotImplementedError

    def calculate_derivative(self, coordinate: List[float], axis: str) -> float:
        """
        Get the derivative of the function at the given coordinate
        """
        raise NotImplementedError

    def calculate_arc_length(self, a: float, b: float) -> float:
        """
        Get the arc length of the function between a and b.
        """
        raise NotImplementedError
