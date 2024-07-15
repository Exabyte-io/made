from typing import List, Optional

import numpy as np
from mat3ra.made.material import Material
from ..third_party import ASEAtoms, ase_make_supercell
from ..utils import decorator_convert_2x2_to_3x3
from ..convert import from_ase, decorator_convert_material_args_kwargs_to_atoms


@decorator_convert_2x2_to_3x3
@decorator_convert_material_args_kwargs_to_atoms
def create_supercell(
    atoms: ASEAtoms, supercell_matrix: Optional[List[List[int]]] = None, scaling_factor: Optional[List[int]] = None
) -> Material:
    """
    Create a supercell of the atoms.

    Args:
        atoms (Material): The atoms to create a supercell of.
        supercell_matrix (List[List[int]]): The supercell matrix (e.g. [[3,0,0],[0,3,0],[0,0,1]]).
        scaling_factor (List[int], optional): The scaling factor instead of matrix (e.g. [3,3,1]).

    Returns:
        Material: The supercell of the atoms.
    """
    if scaling_factor is not None:
        supercell_matrix = np.multiply(scaling_factor, np.eye(3)).tolist()
    supercell_atoms = ase_make_supercell(atoms, supercell_matrix)
    return Material(from_ase(supercell_atoms))
