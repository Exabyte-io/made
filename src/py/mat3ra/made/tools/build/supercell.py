from typing import List, Optional

import numpy as np
from mat3ra.made.material import Material
from ..third_party import ase_make_supercell
from ..utils import decorator_convert_2x2_to_3x3
from ..convert import from_ase, to_ase


@decorator_convert_2x2_to_3x3
def create_supercell(
    material: Material, supercell_matrix: Optional[List[List[int]]] = None, scaling_factor: Optional[List[int]] = None
) -> Material:
    """
    Create a supercell of the atoms.

    Args:
        material (Material): The atoms to create a supercell of.
        supercell_matrix (List[List[int]]): The supercell matrix (e.g. [[3,0,0],[0,3,0],[0,0,1]]).
        scaling_factor (List[int], optional): The scaling factor instead of matrix (e.g. [3,3,1]).

    Returns:
        Material: The supercell of the atoms.
    """
    atoms = to_ase(material)
    if scaling_factor is not None:
        supercell_matrix = np.multiply(scaling_factor, np.eye(3)).tolist()
    supercell_atoms = ase_make_supercell(atoms, supercell_matrix)
    new_material = Material.create(from_ase(supercell_atoms))
    if material.metadata:
        new_material.metadata = material.metadata
    new_material.name = material.name
    return new_material
