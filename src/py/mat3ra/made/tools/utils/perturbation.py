from typing import List, Literal

import numpy as np
from scipy.integrate import quad
from scipy.optimize import root_scalar

from .functions import AXIS_TO_INDEX_MAP, EQUATION_RANGE_COEFFICIENT, FunctionHolder


class PerturbationFunctionHolder(FunctionHolder):
    def calculate_arc_length_equation(self, w_prime: float, w: float) -> float:
        """
        Get the arc length equation for the perturbation function.
        """
        arc_length = quad(
            lambda t: np.sqrt(1 + (self.calculate_derivative([t]) ** 2)),
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

    def calculate_derivative(self, coordinate: List[float]) -> float:
        w = coordinate[AXIS_TO_INDEX_MAP[self.axis]]
        return self.amplitude * 2 * np.pi / self.wavelength * np.cos(2 * np.pi * w / self.wavelength + self.phase)

    def apply_perturbation(self, coordinate: List[float]) -> List[float]:
        return [coordinate[0], coordinate[1], coordinate[2] + self.apply_function(coordinate)]

    def transform_coordinates(self, coordinate: List[float]) -> List[float]:
        index = AXIS_TO_INDEX_MAP[self.axis]

        w = coordinate[index]
        # Find x' such that the integral from 0 to x' equals x
        result = root_scalar(
            self.calculate_arc_length_equation,
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
