from functools import wraps
from typing import Callable, List, Optional

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
from .factories import PerturbationFunctionHolderFactory

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
