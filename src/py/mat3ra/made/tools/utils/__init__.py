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


def decorator_handle_periodic_boundary_conditions(cutoff):
    def decorator(func):
        @wraps(func)
        def wrapper(material, *args, **kwargs):
            # Augment material with periodic images
            augmented_material, original_count = augment_material_with_periodic_images(material, cutoff)

            # Call the original function with augmented material
            result = func(augmented_material, *args, **kwargs)

            # Filter results to include only original atoms
            if isinstance(result, list):  # Assuming result is a list of indices
                result = [idx for idx in result if idx < original_count]
            return result

        return wrapper

    return decorator


def augment_material_with_periodic_images(material, cutoff):
    """Augment the material's dataset by adding atoms from periodic images near boundaries."""
    from ..build.utils import merge_materials
    from ..modify import filter_by_box, translate_by_vector

    material = material.clone()
    original_count = len(material.basis.coordinates.values)
    material_slice_x1 = filter_by_box(material, [0, 0, 0], [cutoff, 1, 1])
    material_slice_x2 = filter_by_box(material, [1 - cutoff, 0, 0], [1, 1, 1])
    translated_material_slice_x1 = translate_by_vector(material_slice_x1, [1, 0, 0])
    translated_material_slice_x2 = translate_by_vector(material_slice_x2, [-1, 0, 0])
    augmented_material_x = merge_materials([material, translated_material_slice_x1, translated_material_slice_x2])

    material_slice_y1 = filter_by_box(augmented_material_x, [-cutoff, 0, 0], [1 + cutoff, cutoff, 1])
    material_slice_y2 = filter_by_box(augmented_material_x, [-cutoff, 1 - cutoff, 0], [1 + cutoff, 1, 1])
    translated_material_slice_y1 = translate_by_vector(material_slice_y1, [0, 1, 0])
    translated_material_slice_y2 = translate_by_vector(material_slice_y2, [0, -1, 0])
    augmented_material_y = merge_materials(
        [augmented_material_x, translated_material_slice_y1, translated_material_slice_y2]
    )

    material_slice_z1 = filter_by_box(augmented_material_y, [-cutoff, -cutoff, 0], [1 + cutoff, 1 + cutoff, cutoff])
    material_slice_z2 = filter_by_box(augmented_material_y, [-cutoff, -cutoff, 1 - cutoff], [1 + cutoff, 1 + cutoff, 1])
    translated_material_slice_z1 = translate_by_vector(material_slice_z1, [0, 0, 1])
    translated_material_slice_z2 = translate_by_vector(material_slice_z2, [0, 0, -1])
    augmented_material = merge_materials(
        [augmented_material_y, translated_material_slice_z1, translated_material_slice_z2]
    )

    return augmented_material, original_count
