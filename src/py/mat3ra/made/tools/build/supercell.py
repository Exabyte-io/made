from typing import List, Optional

import numpy as np
from mat3ra.esse.models.materials_category_components.entities.auxiliary.three_dimensional.supercell_matrix_3d import (
    SupercellMatrix3DSchema,
)

from mat3ra.made.material import Material
from . import BuildMetadata, BaseConfigurationPydantic
from ..operations.core.unary import supercell


class SupercellConfiguration(BaseConfigurationPydantic):
    type: str = "SupercellConfiguration"
    supercell_matrix: SupercellMatrix3DSchema = SupercellMatrix3DSchema(root=np.eye(3).tolist())


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
    BuildMetadata.set_build_metadata(
        new_material,
        configuration=SupercellConfiguration(supercell_matrix=supercell_matrix),
    )
    return new_material
