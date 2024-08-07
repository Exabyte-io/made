from typing import Callable, Dict, List, Literal, Optional, Tuple

from .factories import PerturbationFunctionHolderFactory
from .functions import AXIS_TO_INDEX_MAP


class PerturbationFunctionHolder:
    @staticmethod
    def get_coord_transformation(perturbation_json: dict) -> Callable:
        new_perturbation_json = perturbation_json.copy()
        name = new_perturbation_json.pop("type")
        helper_function = PerturbationFunctionHolderFactory.get_class_by_name(name)
        # TODO: add type of SineWave (or corresponding one to the return of the factory)
        return helper_function.get_transform_coordinates(**new_perturbation_json)

    @staticmethod
    def sine_wave(
        amplitude: float = 0.1,
        wavelength: float = 1,
        phase: float = 0,
        axis: Optional[Literal["x", "y"]] = "x",
    ) -> Tuple[Callable[[List[float]], List[float]], Dict]:
        """
        Deform a coordinate using a sine wave.
        Args:
            amplitude (float): The amplitude of the sine wave in cartesian coordinates.
            wavelength (float): The wavelength of the sine wave in cartesian coordinates.
            phase (float): The phase of the sine wave in cartesian coordinates.
            axis (str): The axis of the direction of the sine wave.

        Returns:
            Tuple[Callable[[List[float]], List[float]], Dict]: The perturbation function and its configuration
        """
        if axis in AXIS_TO_INDEX_MAP:
            index = AXIS_TO_INDEX_MAP[axis]
        perturbation_function = PerturbationFunctionHolderFactory.get_class_by_name("sine_wave")

        def perturbation(coordinate: List[float]):
            return [
                coordinate[0],
                coordinate[1],
                coordinate[2] + perturbation_function.get_function(coordinate[index], amplitude, wavelength, phase),
            ]

        config = perturbation_function.get_json(amplitude, wavelength, phase, axis)
        return perturbation, config
