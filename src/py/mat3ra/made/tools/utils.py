from functools import wraps
from typing import Callable

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
