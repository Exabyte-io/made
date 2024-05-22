from pymatgen.core.structure import Structure
from mat3ra.utils.matrix import convert_2x2_to_3x3


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


def decorator_convert_2x2_to_3x3(func):
    """
    Decorator to convert a 2x2 matrix to a 3x3 matrix.
    Args:
        func: The function to decorate.

    Returns:
        The decorated function.
    """

    def wrapper(matrix):
        return convert_2x2_to_3x3(matrix) if len(matrix) == 2 else matrix

    return wrapper
