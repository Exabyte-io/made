from typing import List, Optional

import numpy as np

from mat3ra.made.material import Material
from ..operations.core.unary import supercell


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
    if scaling_factor is not None:
        supercell_matrix = np.multiply(scaling_factor, np.eye(3)).tolist()
    new_material = supercell(material, supercell_matrix)
    return new_material
