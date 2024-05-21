from typing import TypeVar, List
from mat3ra.utils.matrix import convert_2x2_to_3x3

from ...material import Material
from ..convert import from_ase, to_ase

# Used to allow parameters based on Material (Slab, Interface etc)
MaterialBased = TypeVar("MaterialBased", bound=Material)


def create_supercell(material: Material, supercell_matrix: List[List[int]]) -> Material:
    """
    Create a supercell of the material.

    Args:
        material (Material): The material to create a supercell of.
        supercell_matrix (List[List[int]]): The supercell matrix.

    Returns:
        Material: The supercell of the material.
    """
    from ase.build.supercells import make_supercell

    # if xy_supercell provided, create 3x3 matrix
    if len(supercell_matrix) == 2:
        supercell_matrix = convert_2x2_to_3x3(supercell_matrix)

    supercell_atoms = make_supercell(to_ase(material), supercell_matrix)
    return Material(from_ase(supercell_atoms))
