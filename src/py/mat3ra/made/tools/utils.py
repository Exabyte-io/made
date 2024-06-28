from functools import wraps
from typing import Callable, List

import numpy as np
from mat3ra.made.basis import Basis
from mat3ra.utils.matrix import convert_2x2_to_3x3

from .third_party import PymatgenStructure


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
