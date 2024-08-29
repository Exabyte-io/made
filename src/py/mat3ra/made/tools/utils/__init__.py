from functools import wraps
from typing import Callable, List, Optional

import numpy as np
from mat3ra.made.material import Material
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


def decorator_handle_periodic_boundary_conditions(cutoff):
    def decorator(func):
        @wraps(func)
        def wrapper(material, *args, **kwargs):
            augmented_material, original_count = augment_material_with_periodic_images(material, cutoff)
            result = func(augmented_material, *args, **kwargs)

            if isinstance(result, list):
                result = [idx for idx in result if idx < original_count]

            if isinstance(result, list) and all(isinstance(coord, float) for coord in result):
                result = [coordinate for coordinate in result if 0 <= coordinate <= 1]

            return result

        return wrapper

    return decorator


def filter_and_translate(coordinates: np.ndarray, elements: np.ndarray, axis: int, cutoff: float, direction: int):
    """
    Filter and translate atom coordinates based on the axis and direction.

    Args:
        coordinates (np.ndarray): The coordinates of the atoms.
        elements (np.ndarray): The elements of the atoms.
        axis (int): The axis to filter and translate.
        cutoff (float): The cutoff value for filtering.
        direction (int): The direction to translate.

    Returns:
        Tuple[np.ndarray, np.ndarray]: The filtered and translated coordinates and elements.
    """
    mask = (coordinates[:, axis] < cutoff) if direction == 1 else (coordinates[:, axis] > (1 - cutoff))
    filtered_coordinates = coordinates[mask]
    filtered_elements = elements[mask]
    translation_vector = np.zeros(3)
    translation_vector[axis] = direction
    translated_coordinates = filtered_coordinates + translation_vector
    return translated_coordinates, filtered_elements


def augment_material_with_periodic_images(material: Material, cutoff: float = 0.1):
    """
    Augment the material's dataset by adding atoms from periodic images near boundaries.

    Args:
        material (Material): The material to augment.
        cutoff (float): The cutoff value for filtering atoms near boundaries.

    Returns:
        Tuple[Material, int]: The augmented material and the original count of atoms.
    """
    original_count = len(material.basis.coordinates.values)
    coordinates = np.array(material.basis.coordinates.values)
    elements = np.array(material.basis.elements.values)
    augmented_coords = coordinates.tolist()
    augmented_elems = elements.tolist()

    for axis in range(3):
        for direction in [-1, 1]:
            translated_coords, translated_elems = filter_and_translate(coordinates, elements, axis, cutoff, direction)
            augmented_coords.extend(translated_coords)
            augmented_elems.extend(translated_elems)
            coordinates = np.array(augmented_coords)
            elements = np.array(augmented_elems)

    augmented_material = material.clone()
    new_basis = augmented_material.basis.copy()
    for i, coord in enumerate(augmented_coords):
        new_basis.add_atom(augmented_elems[i], coord)

    augmented_material.basis = new_basis
    return augmented_material, original_count
