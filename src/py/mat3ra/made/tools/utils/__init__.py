from functools import wraps
from typing import Any, Callable, List, Optional, Union, cast

import numpy as np
from mat3ra.esse.models.material.reusable.supercell_matrix_2d import SupercellMatrix2DSchema
from mat3ra.esse.models.materials_category.single_material.two_dimensional.slab.configuration import (
    SupercellMatrix2DSchemaItem,
)
from mat3ra.utils.matrix import convert_2x2_to_3x3

DEFAULT_SCALING_FACTOR = np.array([3, 3, 3])
DEFAULT_TRANSLATION_VECTOR = 1 / DEFAULT_SCALING_FACTOR


def is_primitive_2x2_matrix(matrix: Any) -> bool:
    return (
        isinstance(matrix, list) and len(matrix) == 2 and all(isinstance(row, list) and len(row) == 2 for row in matrix)
    )


def normalize_2x2_matrix(
    matrix: Union[
        List[List[float]],
        SupercellMatrix2DSchema,
    ]
) -> Optional[List[List[float]]]:
    """
    Normalize any matrix-like structure to a plain 2x2 list of floats.
    Returns None if normalization is not possible.
    """
    if isinstance(matrix, SupercellMatrix2DSchema) and matrix.root:
        matrix = matrix.root

    if isinstance(matrix, SupercellMatrix2DSchemaItem):
        return [matrix.root]

    if (
        isinstance(matrix, list)
        and len(matrix) == 2
        and all(isinstance(row, SupercellMatrix2DSchemaItem) for row in matrix)
    ):
        return [cast(SupercellMatrix2DSchemaItem, row).root for row in matrix]

    if is_primitive_2x2_matrix(matrix):
        return matrix  # already normalized

    return None  # unrecognized format


def decorator_convert_2x2_to_3x3(func: Callable) -> Callable:
    """
    Decorator that converts a 2x2 matrix input to a 3x3 matrix.
    Supports schema-based formats and raw nested lists.
    """

    @wraps(func)
    def wrapper(*args, **kwargs):
        def convert_if_matrix(arg):
            matrix = normalize_2x2_matrix(arg)
            return convert_2x2_to_3x3(matrix) if matrix else arg

        args = tuple(convert_if_matrix(arg) for arg in args)

        if "supercell_matrix" in kwargs:
            kwargs["supercell_matrix"] = convert_if_matrix(kwargs["supercell_matrix"])

        return func(*args, **kwargs)

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
