from functools import wraps
from typing import Callable, List, Optional

import numpy as np
from mat3ra.made.material import Material
from mat3ra.made.utils import ArrayWithIds
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
    """
    Decorator to handle periodic boundary conditions.

    Copies atoms near boundaries within the cutoff distance beyond the opposite side of the cell
    creating the effect of periodic boundary conditions for edge atoms.

    Results of the function are filtered to remove atoms or coordinates outside the original cell.

    Args:
        cutoff (float): The cutoff distance for a border slice in crystal coordinates.

    Returns:
        Callable: The decorated function.
    """

    def decorator(func):
        @wraps(func)
        def wrapper(material, *args, **kwargs):
            augmented_material, last_id = augment_material_with_periodic_images(material, cutoff)
            result = func(augmented_material, *args, **kwargs)

            if isinstance(result, list):
                if all(isinstance(x, int) for x in result):
                    result = [id for id in result if id <= last_id]
                elif all(isinstance(x, list) and len(x) == 3 for x in result):
                    result = [coord for coord in result if all(0 <= c < 1 for c in coord)]
            return result

        return wrapper

    return decorator


def augment_material_with_periodic_images(material: Material, cutoff: float = 0.1):
    """
    Augment the material's dataset by adding atoms from periodic images within a cutoff distance from the boundaries by
    copying them to the opposite side of the cell, translated by the cell vector beyond the boundary.

    Args:
        material (Material): The material to augment.
        cutoff (float): The cutoff value for filtering atoms near boundaries, in crystal coordinates.

    Returns:
        Tuple[Material, int]: The augmented material and the original count of atoms.
    """
    last_id = material.basis.coordinates.ids[-1]
    coordinates = np.array(material.basis.coordinates.values)
    elements = np.array(material.basis.elements.values)
    augmented_material = material.clone()
    new_basis = augmented_material.basis.copy()

    for axis in range(3):
        for direction in [-1, 1]:
            mask = (coordinates[:, axis] < cutoff) if direction == 1 else (coordinates[:, axis] > (1 - cutoff))
            filtered_coordinates = coordinates[mask]
            filtered_elements = elements[mask]
            translation_vector = np.zeros(3)
            translation_vector[axis] = direction
            translated_coordinates = filtered_coordinates + translation_vector
            for coord, elem in zip(translated_coordinates, filtered_elements):
                new_basis.add_atom(elem, coord)

    augmented_material.basis = new_basis
    return augmented_material, last_id
