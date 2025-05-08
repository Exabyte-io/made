from typing import Optional

import numpy as np
from mat3ra.esse.models.materials_category.single_material.two_dimensional.slab.configuration import (
    SupercellMatrix2DSchemaItem,
)
from mat3ra.utils.matrix import convert_2x2_to_3x3 as convert_2x2_to_3x3_orig

DEFAULT_SCALING_FACTOR = np.array([3, 3, 3])
DEFAULT_TRANSLATION_VECTOR = 1 / DEFAULT_SCALING_FACTOR


from functools import wraps
from typing import Callable, List


def convert_2x2_to_3x3(matrix):
    """
    Convert a 2x2 matrix (as list of lists, list of SupercellMatrix2DSchemaItem, or SupercellMatrix2DSchemaItem)
    to a 3x3 matrix by adding a third unitary orthogonal basis vector.
    """
    # If it's a list of SupercellMatrix2DSchemaItem, extract .root from each
    if (
        isinstance(matrix, list)
        and len(matrix) == 2
        and all(isinstance(row, SupercellMatrix2DSchemaItem) for row in matrix)
    ):
        matrix = [row.root for row in matrix]
    elif isinstance(matrix, SupercellMatrix2DSchemaItem):
        matrix = matrix.root
    return convert_2x2_to_3x3_orig(matrix)


def decorator_convert_2x2_to_3x3(func: Callable) -> Callable:
    """
    Decorator that converts a 2x2 matrix input to a 3x3 matrix.
    Supports:
    - SupercellMatrix2DSchema
    - List[List[float]] with 2x2 shape
    """

    @wraps(func)
    def wrapper(*args, **kwargs):
        def needs_conversion(arg):
            # list of SupercellMatrix2DSchemaItem
            if (
                isinstance(arg, list)
                and len(arg) == 2
                and all(isinstance(row, SupercellMatrix2DSchemaItem) for row in arg)
            ):
                return True
            # single SupercellMatrix2DSchemaItem
            if isinstance(arg, SupercellMatrix2DSchemaItem):
                return True
            # plain 2x2 list
            if isinstance(arg, list) and len(arg) == 2 and all(isinstance(row, list) and len(row) == 2 for row in arg):
                return True
            return False

        new_args = [convert_2x2_to_3x3(arg) if needs_conversion(arg) else arg for arg in args]
        return func(*new_args, **kwargs)

    return wrapper


def decorator_convert_position_to_coordinate(func: Callable) -> Callable:
    """
    A decorator that converts a 2D position [x, y] to a 3D coordinate [x, y, 0.0].
    """

    @wraps(func)
    def wrapper(*args, **kwargs):
        if "coordinate" in kwargs:
            coordinate = kwargs["coordinate"]
        else:
            # Find the position of the 'coordinate' argument and get its value
            coordinate_index = func.__code__.co_varnames.index("coordinate")
            coordinate = args[coordinate_index]

        if len(coordinate) == 2:
            coordinate = [coordinate[0], coordinate[1], 0.0]

        if "coordinate" in kwargs:
            kwargs["coordinate"] = coordinate
        else:
            args = list(args)
            args[coordinate_index] = coordinate
            args = tuple(args)

        return func(*args, **kwargs)

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
