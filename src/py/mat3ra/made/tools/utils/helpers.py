from typing import Callable, List, Literal

import numpy as np
from pydantic import BaseModel
from scipy.integrate import quad
from scipy.optimize import root_scalar

AXIS_TO_INDEX_MAP = {"x": 0, "y": 1, "z": 2}
EQUATION_RANGE_COEFFICIENT = 5


class FunctionHolder(BaseModel):

    @staticmethod
    def get_function(*args, **kwargs):
        """
        Get the function of the perturbation.
        """
        raise NotImplementedError

    @staticmethod
    def get_derivative(*args, **kwargs):
        """
        Get the derivative of the perturbation function.
        """
        raise NotImplementedError

    @staticmethod
    def get_arc_length_equation(*args, **kwargs):
        """
        Get the equation to calculate the arc length between [0,0,0] and a given coordinate of the perturbation.
        """
        raise NotImplementedError

    @staticmethod
    def get_transform_coordinates(*args, **kwargs):
        """
        Transform coordinates to preserve the distance between points on a sine wave when perturbation is applied.
        Achieved by calculating the integral of the length between [0,0,0] and given coordinate.

        Returns:
            Callable[[List[float]], List[float]]: The coordinates transformation function.
        """
        raise NotImplementedError


class SineWave(FunctionHolder):

    @staticmethod
    def get_function(
        w: float, amplitude: float, wavelength: float, phase: float
    ) -> Callable[[List[float]], List[float]]:
        return amplitude * np.sin(2 * np.pi * w / wavelength + phase)

    @staticmethod
    def get_derivative(w: float, amplitude: float, wavelength: float, phase: float) -> float:
        return amplitude * 2 * np.pi / wavelength * np.cos(2 * np.pi * w / wavelength + phase)

    @staticmethod
    def get_arc_length_equation(w_prime: float, w: float, amplitude: float, wavelength: float, phase: float) -> float:
        arc_length = quad(
            lambda t: np.sqrt(1 + (SineWave.get_derivative(t, amplitude, wavelength, phase)) ** 2),
            a=0,
            b=w_prime,
        )[0]
        return arc_length - w

    @staticmethod
    def get_transform_coordinates(
        amplitude: float, wavelength: float, phase: float, axis: Literal["x", "y"]
    ) -> Callable[[List[float]], List[float]]:
        index = AXIS_TO_INDEX_MAP[axis]

        def coordinate_transformation(coordinate: List[float]):
            w = coordinate[index]
            # Find x' such that the integral from 0 to x' equals x
            result = root_scalar(
                SineWave.get_arc_length_equation,
                args=(w, amplitude, wavelength, phase),
                bracket=[0, EQUATION_RANGE_COEFFICIENT * w],
                method="brentq",
            )
            coordinate[index] = result.root
            return coordinate

        return coordinate_transformation
