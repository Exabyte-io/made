from typing import Callable, List, Literal

import numpy as np
from pydantic import BaseModel
from scipy.integrate import quad
from scipy.optimize import root_scalar

AXIS_TO_INDEX_MAP = {"x": 0, "y": 1, "z": 2}
EQUATION_RANGE_COEFFICIENT = 5


class FunctionHolder(BaseModel):
    @staticmethod
    def get_derivative(*args, **kwargs):
        raise NotImplementedError

    @staticmethod
    def get_arc_length_equation(*args, **kwargs):
        raise NotImplementedError

    @staticmethod
    def get_transform_coordinates(*args, **kwargs):
        raise NotImplementedError


class SineWave(FunctionHolder):
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
        """
        Transform coordinates to preserve the distance between points on a sine wave.
        Achieved by calculating the integral of the length between [0,0,0] and given coordinate.

        Args:
            amplitude (float): The amplitude of the sine wave in cartesian coordinates.
            wavelength (float): The wavelength of the sine wave in cartesian coordinates.
            phase (float): The phase of the sine wave in cartesian coordinates.
            axis (str): The axis of the direction of the sine wave.

        Returns:
            Callable[[List[float]], List[float]]: The coordinates transformation function.
        """
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
