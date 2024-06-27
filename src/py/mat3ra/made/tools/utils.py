from functools import wraps
from typing import Callable, List

import numpy as np
from mat3ra.made.basis import Basis
from mat3ra.utils.matrix import convert_2x2_to_3x3
from pymatgen.core.structure import Structure


# TODO: convert to accept ASE Atoms object
def translate_to_bottom_pymatgen_structure(structure: Structure):
    """
    Translate the structure to the bottom of the cell.
    Args:
        structure (Structure): The pymatgen Structure object to translate.

    Returns:
        Structure: The translated pymatgen Structure object.
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


def convert_basis_to_cartesian(basis: Basis) -> Basis:
    """
    Convert the basis to the Cartesian coordinates.
    Args:
        basis (Dict): The basis to convert.

    Returns:
        Dict: The basis in Cartesian coordinates.
    """
    if basis.units == "cartesian":
        return basis
    unit_cell = np.array(basis.cell)
    basis.coordinates = np.multiply(basis.coordinates, unit_cell)
    basis.units = "cartesian"
    return basis


def convert_basis_to_crystal(basis: Basis) -> Basis:
    """
    Convert the basis to the crystal coordinates.
    Args:
        basis (Dict): The basis to convert.

    Returns:
        Dict: The basis in crystal coordinates.
    """
    if basis.units == "crystal":
        return basis
    unit_cell = np.array(basis.cell)
    basis.coordinates.values = np.multiply(basis.coordinates.values, np.linalg.inv(unit_cell))
    basis.units = "crystal"
    return basis


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
# TODO: Predefined equations can be exported using a factory or enum
def is_point_in_circle(x=0, y=0, r=1) -> Callable[[List[float]], bool]:
    """
    Check if a point is inside a circle.
    Args:
        x (float): The x-coordinate of the circle center.
        y (float): The y-coordinate of the circle center.
        r (float): The radius of the circle.

    Returns:
        Callable[[List[float]], bool]: The condition function to check if a point is inside the circle.
    """
    return lambda vector: (vector[0] - x) ** 2 + (vector[1] - y) ** 2 <= r**2


def is_point_in_rectangle(a=0, b=1, c=0, d=1) -> Callable[[List[float]], bool]:
    """
    Check if a point is inside a rectangle.
    Args:
        a (float): Lower limit of x-coordinate.
        b (float): Upper limit of x-coordinate.
        c (float): Lower limit of y-coordinate.
        d (float): Upper limit of y-coordinate.

    Returns:
        Callable[[List[float]], bool]: The condition function to check if a point is inside the rectangle.
    """
    return lambda vector: a <= vector[0] <= b and c <= vector[1] <= d
