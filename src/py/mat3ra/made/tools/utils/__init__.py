from functools import wraps
from typing import Callable, Dict, List, Literal, Optional, Tuple

import numpy as np
from mat3ra.utils.matrix import convert_2x2_to_3x3

from ..third_party import PymatgenStructure
from .coordinate import (
    is_coordinate_behind_plane,
    is_coordinate_in_box,
    is_coordinate_in_cylinder,
    is_coordinate_in_sphere,
    is_coordinate_in_triangular_prism,
)
from .factories import PerturbationFunctionHelperFactory
from .helpers import AXIS_TO_INDEX_MAP

DEFAULT_SCALING_FACTOR = np.array([3, 3, 3])
DEFAULT_TRANSLATION_VECTOR = 1 / DEFAULT_SCALING_FACTOR


# TODO: convert to accept ASE Atoms object
def translate_to_bottom_pymatgen_structure(structure: PymatgenStructure):
    """
    Translate the structure to the bottom of the cell.
    Args:
        structure (PymatgenStructure): The pymatgen Structure object to translate.

    Returns:
        PymatgenStructure: The translated pymatgen Structure object.
    """
    min_c = min(site.c for site in structure)
    translation_vector = [0, 0, -min_c]
    translated_structure = structure.copy()
    for site in translated_structure:
        site.coords += translation_vector
    return translated_structure


def decorator_convert_2x2_to_3x3(func: Callable) -> Callable:
    """
    Decorator to convert a 2x2 matrix to a 3x3 matrix.
    """

    @wraps(func)
    def wrapper(*args, **kwargs):
        new_args = [convert_2x2_to_3x3(arg) if isinstance(arg, list) and len(arg) == 2 else arg for arg in args]
        return func(*new_args, **kwargs)

    return wrapper


def get_distance_between_coordinates(coordinate1: List[float], coordinate2: List[float]) -> float:
    """
    Get the distance between two coordinates.
    Args:
        coordinate1 (List[float]): The first coordinate.
        coordinate2 (List[float]): The second coordinate.

    Returns:
        float: The distance between the two coordinates.
    """
    return float(np.linalg.norm(np.array(coordinate1) - np.array(coordinate2)))


def get_norm(vector: List[float]) -> float:
    """
    Get the norm of a vector.
    Args:
        vector (List[float]): The vector.

    Returns:
        float: The norm of the vector.
    """
    return float(np.linalg.norm(vector))


# Condition functions:


def transform_coordinate_to_supercell(
    coordinate: List[float],
    scaling_factor: Optional[List[int]] = None,
    translation_vector: Optional[List[float]] = None,
    reverse: bool = False,
) -> List[float]:
    """
    Convert a crystal coordinate of unit cell to a coordinate in a supercell.
    Args:
        coordinate (List[float]): The coordinates to convert.
        scaling_factor (List[int]): The scaling factor for the supercell.
        translation_vector (List[float]): The translation vector for the supercell.
        reverse (bool): Whether to convert in the reverse transformation.

    Returns:
        List[float]: The converted coordinates.
    """
    if scaling_factor is None:
        np_scaling_factor = np.array([3, 3, 3])
    else:
        np_scaling_factor = np.array(scaling_factor)

    if translation_vector is None:
        np_translation_vector = np.array([0, 0, 0])
    else:
        np_translation_vector = np.array(translation_vector)

    np_coordinate = np.array(coordinate)
    converted_array = np_coordinate * (1 / np_scaling_factor) + np_translation_vector
    if reverse:
        converted_array = (np_coordinate - np_translation_vector) * np_scaling_factor
    return converted_array.tolist()


class CoordinateConditionBuilder:
    @staticmethod
    def create_condition(condition_type: str, evaluation_func: Callable, **kwargs) -> Tuple[Callable, Dict]:
        condition_json = {"type": condition_type, **kwargs}
        return lambda coordinate: evaluation_func(coordinate, **kwargs), condition_json

    @staticmethod
    def cylinder(center_position=None, radius: float = 0.25, min_z: float = 0, max_z: float = 1):
        if center_position is None:
            center_position = [0.5, 0.5]
        return CoordinateConditionBuilder.create_condition(
            condition_type="cylinder",
            evaluation_func=is_coordinate_in_cylinder,
            center_position=center_position,
            radius=radius,
            min_z=min_z,
            max_z=max_z,
        )

    @staticmethod
    def sphere(center_position=None, radius: float = 0.25):
        if center_position is None:
            center_position = [0.5, 0.5, 0.5]
        return CoordinateConditionBuilder.create_condition(
            condition_type="sphere",
            evaluation_func=is_coordinate_in_sphere,
            center_position=center_position,
            radius=radius,
        )

    @staticmethod
    def triangular_prism(
        position_on_surface_1: List[float] = [0, 0],
        position_on_surface_2: List[float] = [1, 0],
        position_on_surface_3: List[float] = [0, 1],
        min_z: float = 0,
        max_z: float = 1,
    ):
        return CoordinateConditionBuilder.create_condition(
            condition_type="prism",
            evaluation_func=is_coordinate_in_triangular_prism,
            coordinate_1=position_on_surface_1,
            coordinate_2=position_on_surface_2,
            coordinate_3=position_on_surface_3,
            min_z=min_z,
            max_z=max_z,
        )

    @staticmethod
    def box(min_coordinate=None, max_coordinate=None):
        if max_coordinate is None:
            max_coordinate = [1, 1, 1]
        if min_coordinate is None:
            min_coordinate = [0, 0, 0]
        return CoordinateConditionBuilder.create_condition(
            condition_type="box",
            evaluation_func=is_coordinate_in_box,
            min_coordinate=min_coordinate,
            max_coordinate=max_coordinate,
        )

    @staticmethod
    def plane(plane_normal: List[float], plane_point_coordinate: List[float]):
        return CoordinateConditionBuilder.create_condition(
            condition_type="plane",
            evaluation_func=is_coordinate_behind_plane,
            plane_normal=plane_normal,
            plane_point_coordinate=plane_point_coordinate,
        )


class PerturbationFunctionHolder:
    @staticmethod
    def get_coord_transformation(perturbation_json: dict) -> Callable:
        new_perturbation_json = perturbation_json.copy()
        name = new_perturbation_json.pop("type")
        helper_function = PerturbationFunctionHelperFactory.get_class_by_name(name)
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
        perturbation_function = PerturbationFunctionHelperFactory.get_class_by_name("sine_wave")

        def perturbation(coordinate: List[float]):
            return [
                coordinate[0],
                coordinate[1],
                coordinate[2] + perturbation_function.get_function(coordinate[index], amplitude, wavelength, phase),
            ]

        config = {"type": "sine_wave", "amplitude": amplitude, "wavelength": wavelength, "phase": phase, "axis": axis}

        return perturbation, config

    @staticmethod
    def sine_wave_radial(
        amplitude: float = 0.1, wavelength: float = 1, phase: float = 0, center_position=None
    ) -> Tuple[Callable[[List[float]], List[float]], Dict]:
        """
        Deform a coordinate using a radial sine wave.
        Args:
            amplitude (float): The amplitude of the sine wave in cartesian coordinates.
            wavelength (float): The wavelength of the sine wave in cartesian coordinates.
            phase (float): The phase of the sine wave in cartesian coordinates.
            center_position (List[float]): The center position of the sine wave on the plane.

        Returns:
            Tuple[Callable[[List[float]], List[float]], Dict]: The perturbation function and its configuration
        """
        if center_position is None:
            center_position = [0.5, 0.5]

        def perturbation(coordinate: List[float]):
            np_position = np.array(coordinate[:2])
            np_center_position = np.array(center_position)
            distance = np.linalg.norm(np_position - np_center_position)
            return [
                coordinate[0],
                coordinate[1],
                coordinate[2] + amplitude * np.sin(2 * np.pi * distance / wavelength + phase),
            ]

        config = {
            "type": "sine_wave_radial",
            "amplitude": amplitude,
            "wavelength": wavelength,
            "phase": phase,
            "center_position": center_position,
        }

        return perturbation, config
