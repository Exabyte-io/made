from typing import Any, Callable, List, Optional

import numpy as np
import sympy as sp
from scipy.integrate import quad
from scipy.optimize import root_scalar

from .functions import AXIS_TO_INDEX_MAP, EQUATION_RANGE_COEFFICIENT, FunctionHolder


def default_function() -> sp.Expr:
    return sp.Symbol("f")


class PerturbationFunctionHolder(FunctionHolder):
    variables: List[str] = ["x"]
    symbols: List[sp.Symbol] = [sp.Symbol(var) for var in variables]
    function: sp.Expr = sp.Symbol("f")
    function_numeric: Callable = default_function
    derivatives_numeric: dict = {}

    class Config:
        arbitrary_types_allowed = True

    def __init__(self, function: Optional[sp.Expr] = None, variables: Optional[List[str]] = None, **data: Any):
        """
        Initializes with a function involving multiple variables.
        """
        if function is None:
            function = default_function()
        if variables is None:
            variables = ["x", "y", "z"]
        super().__init__(**data)
        self.variables = variables
        self.symbols = sp.symbols(variables)
        self.function = function
        self.function_numeric = sp.lambdify(self.symbols, self.function, modules=["numpy"])
        self.derivatives_numeric = {
            var: sp.lambdify(self.symbols, sp.diff(self.function, var), modules=["numpy"]) for var in variables
        }

    def apply_function(self, coordinate: List[float]) -> float:
        values = [coordinate[AXIS_TO_INDEX_MAP[var]] for var in self.variables]
        return self.function_numeric(*values)

    def apply_derivative(self, coordinate: List[float], axis: str) -> float:
        if axis in self.variables:
            values = [coordinate[AXIS_TO_INDEX_MAP[var]] for var in self.variables]
            return self.derivatives_numeric[axis](*values)
        else:
            return 0

    def get_arc_length_equation(self, w_prime: float, coordinate: List[float], axis: str) -> float:
        """
        Calculate arc length considering a change along one specific axis.
        """
        index = AXIS_TO_INDEX_MAP[axis]
        a, b = 0, w_prime

        def integrand(t):
            temp_coordinate = coordinate[:]
            temp_coordinate[index] = t
            return np.sqrt(1 + self.apply_derivative(temp_coordinate, axis) ** 2)

        arc_length = quad(integrand, a, b)[0]
        return arc_length - coordinate[index]

    def transform_coordinates(self, coordinate: List[float]) -> List[float]:
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

    def apply_perturbation(self, coordinate: List[float]) -> List[float]:
        """
        Apply the perturbation to the given coordinate by adding the function's value to the third coordinate (z-axis).
        """
        perturbation_value = self.apply_function(coordinate)
        perturbed_coordinate = coordinate[:]
        perturbed_coordinate[2] += perturbation_value
        return perturbed_coordinate

    def get_json(self) -> dict:
        return {
            "type": self.__class__.__name__,
            "function": str(self.function),
            "variables": self.variables,
        }


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
        function = self._create_function(amplitude, wavelength, phase, axis)
        variables = [axis]
        super().__init__(function=function, variables=variables, **data)
        self.amplitude = amplitude
        self.wavelength = wavelength
        self.phase = phase
        self.axis = axis

    def _create_function(self, amplitude, wavelength, phase, axis) -> sp.Expr:
        w = sp.Symbol(axis)
        return amplitude * sp.sin(2 * sp.pi * w / wavelength + phase)

    def get_json(self) -> dict:
        return {
            "type": self.__class__.__name__,
            "function": str(self.function),
            "variables": self.variables,
            "amplitude": self.amplitude,
            "wavelength": self.wavelength,
            "phase": self.phase,
            "axis": self.axis,
        }
