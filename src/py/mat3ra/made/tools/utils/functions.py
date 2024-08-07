from typing import List, Literal

import numpy as np
from pydantic import BaseModel
from scipy.integrate import quad
from scipy.optimize import root_scalar

AXIS_TO_INDEX_MAP = {"x": 0, "y": 1, "z": 2}
EQUATION_RANGE_COEFFICIENT = 5


class FunctionHolder(BaseModel):
    def apply_function(self, coordinate: List[float]) -> float:
        """
        Get the value of the function at the given coordinate.
        """
        raise NotImplementedError

    def apply_derivative(self, coordinate: List[float]) -> float:
        """
        Get the derivative of the function at the given coordinate
        """
        raise NotImplementedError

    def get_arc_length(self, a: float, b: float) -> float:
        """
        Get the arc length of the function between a and b.
        """
        raise NotImplementedError

    def get_json(self) -> dict:
        """
        Get the json representation of the function holder.
        """
        raise NotImplementedError


class PerturbationFunctionHolder(FunctionHolder):
    def get_arc_length_equation(self, w_prime: float, w: float) -> float:
        """
        Get the arc length equation for the perturbation function.
        """
        arc_length = quad(
            lambda t: np.sqrt(1 + (self.apply_derivative([t]) ** 2)),
            a=0,
            b=w_prime,
        )[0]
        return arc_length - w

    def transform_coordinates(self, coordinate: List[float]) -> List[float]:
        """
        Transform coordinates to preserve the distance between points on a sine wave when perturbation is applied.
        Achieved by calculating the integral of the length between [0,0,0] and given coordinate.

        Returns:
            Callable[[List[float]], List[float]]: The coordinates transformation function.
        """
        raise NotImplementedError

    def apply_perturbation(self, coordinate: List[float]) -> List[float]:
        """
        Apply the perturbation to the given coordinate.
        """
        raise NotImplementedError


class SineWavePerturbationFunctionHolder(PerturbationFunctionHolder):
    amplitude: float = 0.05
    wavelength: float = 1
    phase: float = 0
    axis: Literal["x", "y"] = "x"

    def apply_function(self, coordinate: List[float]) -> float:
        w = coordinate[AXIS_TO_INDEX_MAP[self.axis]]
        return self.amplitude * np.sin(2 * np.pi * w / self.wavelength + self.phase)

    def apply_derivative(self, coordinate: List[float]) -> float:
        w = coordinate[AXIS_TO_INDEX_MAP[self.axis]]
        return self.amplitude * 2 * np.pi / self.wavelength * np.cos(2 * np.pi * w / self.wavelength + self.phase)

    def apply_perturbation(self, coordinate: List[float]) -> List[float]:
        return [coordinate[0], coordinate[1], coordinate[2] + self.apply_function(coordinate)]

    def transform_coordinates(self, coordinate: List[float]) -> List[float]:
        index = AXIS_TO_INDEX_MAP[self.axis]

        w = coordinate[index]
        # Find x' such that the integral from 0 to x' equals x
        result = root_scalar(
            self.get_arc_length_equation,
            args=w,
            bracket=[0, EQUATION_RANGE_COEFFICIENT * w],
            method="brentq",
        )
        coordinate[index] = result.root
        return coordinate

    def get_json(self) -> dict:
        return {
            "type": self.__class__.__name__,
            "amplitude": self.amplitude,
            "wavelength": self.wavelength,
            "phase": self.phase,
            "axis": self.axis,
        }
