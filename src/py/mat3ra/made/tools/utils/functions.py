from typing import Any, Callable, List, Optional

import numpy as np
import sympy as sp
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

    def apply_derivative(self, coordinate: List[float], axis: str) -> float:
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


def default_function(coordinate: List[float]) -> float:
    return 0


class PerturbationFunctionHolder(FunctionHolder):
    variables: List[str] = ["x"]
    symbols: List[sp.Symbol] = [sp.Symbol(var) for var in variables]
    function: sp.Expr = sp.Symbol("f")
    function_numeric: Callable = default_function
    derivatives_numeric: dict = {}

    class Config:
        arbitrary_types_allowed = True

    def __init__(self, function: Optional[Callable] = None, variables: Optional[List[str]] = None, **data: Any):
        """
        Initializes with a function involving multiple variables.
        """
        if function is None:
            function = default_function
        if variables is None:
            variables = ["x"]
        super().__init__(**data)
        self.variables = variables
        self.symbols = sp.symbols(variables)
        self.function = function(*self.symbols)
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
        a, b = 0, w_prime  # Integration limits based on the current position along the axis

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
        super().__init__(**data)
        self.amplitude = amplitude
        self.wavelength = wavelength
        self.phase = phase
        self.axis = axis
        function = lambda x: self.amplitude * sp.sin(2 * sp.pi * x / self.wavelength + self.phase)
        variables = [self.axis]

        PerturbationFunctionHolder.__init__(self, function=function, variables=variables)

    def get_json(self) -> dict:
        return {
            "type": self.__class__.__name__,
            "amplitude": self.amplitude,
            "wavelength": self.wavelength,
            "phase": self.phase,
            "axis": self.axis,
        }
