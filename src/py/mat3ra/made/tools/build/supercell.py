from typing import List
from ase import Atoms
from ase.build.supercells import make_supercell

from mat3ra.made.material import Material
from ..utils import decorator_convert_2x2_to_3x3
from ..convert import from_ase, decorator_convert_material_args_kwargs_to_atoms


@decorator_convert_2x2_to_3x3
@decorator_convert_material_args_kwargs_to_atoms
def create_supercell(atoms: Atoms, supercell_matrix: List[List[int]]) -> Material:
    """
    Create a supercell of the atoms.

    Args:
        atoms (Material): The atoms to create a supercell of.
        supercell_matrix (List[List[int]]): The supercell matrix.

    Returns:
        Material: The supercell of the atoms.
    """

    supercell_atoms = make_supercell(atoms, supercell_matrix)
    return Material(from_ase(supercell_atoms))
