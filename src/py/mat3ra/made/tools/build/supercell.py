from typing import List, Optional, Union

import numpy as np
from mat3ra.esse.models.materials_category_components.entities.auxiliary.three_dimensional.supercell_matrix_3d import (
    SupercellMatrix3DSchema,
)

from mat3ra.made.material import Material
from mat3ra.made.tools.build import MaterialWithBuildMetadata, BaseConfigurationPydantic
from ..operations.core.unary import supercell


class SupercellConfiguration(BaseConfigurationPydantic):
    type: str = "SupercellConfiguration"
    supercell_matrix: SupercellMatrix3DSchema = SupercellMatrix3DSchema(root=np.eye(3).tolist())


def create_supercell(
    material: Union[Material, MaterialWithBuildMetadata],
    supercell_matrix: Optional[List[List[int]]] = None,
    scaling_factor: Optional[List[int]] = None,
) -> Material:
    """
    Create a supercell of the atoms.

    Args:
        material (Material): The atoms to create a supercell of.
        supercell_matrix (List[List[int]]): The supercell matrix (e.g. [[3,0,0],[0,3,0],[0,0,1]]).
        scaling_factor (List[int], optional): The scaling factor instead of matrix (e.g. [3,3,1]).
            If provided, supercell_matrix is ignored.

    Returns:
        Material: The supercell material.
    """
    if scaling_factor is not None:
        supercell_matrix = np.multiply(scaling_factor, np.eye(3)).tolist()
    new_material = supercell(material, supercell_matrix)

    # Convert to MaterialWithBuildMetadata and set metadata
    material_with_metadata = MaterialWithBuildMetadata.create(new_material.to_dict())
    configuration = SupercellConfiguration(supercell_matrix=supercell_matrix)
    material_with_metadata.metadata.add_build_metadata_step(
        configuration=configuration.to_dict(), build_parameters={"supercell_matrix": supercell_matrix}
    )
    return material_with_metadata
