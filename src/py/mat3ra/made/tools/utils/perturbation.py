from typing import Any, List

import numpy as np
import sympy as sp
from scipy.integrate import quad
from scipy.optimize import root_scalar

from ...utils import AXIS_TO_INDEX_MAP
from .functions import FunctionHolder

EQUATION_RANGE_COEFFICIENT = 5


class PerturbationFunctionHolder(FunctionHolder):
    def _integrand(self, t: float, coordinate: List[float], axis: str) -> float:
        temp_coordinate = coordinate[:]
        temp_coordinate[AXIS_TO_INDEX_MAP[axis]] = t
        return np.sqrt(1 + self.calculate_derivative(temp_coordinate, axis) ** 2)

    def get_arc_length_equation(self, w_prime: float, coordinate: List[float], axis: str) -> float:
        """
        Calculate arc length considering a change along one specific axis.
        """
        a, b = 0, w_prime
        arc_length = quad(self._integrand, a, b, args=(coordinate, axis), limit=1000)[0]
        return arc_length - coordinate[AXIS_TO_INDEX_MAP[axis]]

    def normalize_coordinates(self, coordinate: List[float]) -> List[float]:
        """
        Transform coordinates to preserve the distance between points on a sine wave when perturbation is applied.
        Achieved by calculating the integral of the length between [0,0,0] and given coordinate.

        Returns:
            Callable[[List[float]], List[float]]: The coordinates transformation function.
        """
        for i, var in enumerate(self.variables):
            index = AXIS_TO_INDEX_MAP[var]
            w = coordinate[index]
            result = root_scalar(
                self.get_arc_length_equation,
                args=(coordinate, var),
                bracket=[0, EQUATION_RANGE_COEFFICIENT * w],
                method="brentq",
            )
            coordinate[index] = result.root
        return coordinate


class SineWavePerturbationFunctionHolder(PerturbationFunctionHolder):
    amplitude: float = 0.05
    wavelength: float = 1
    phase: float = 0
    axis: str = "x"

    def __init__(
        self,
        amplitude: float = 0.05,
        wavelength: float = 1,
        phase: float = 0,
        axis: str = "x",
        **data: Any,
    ):
        w = sp.Symbol(axis)
        function = amplitude * sp.sin(2 * sp.pi * w / wavelength + phase)
        variables = [axis]
        super().__init__(function=function, variables=variables, **data)
        self.amplitude = amplitude
        self.wavelength = wavelength
        self.phase = phase
        self.axis = axis
