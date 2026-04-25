from typing import Optional, Union

from mat3ra.esse.models.core.abstract.matrix_3x3 import Matrix3x3Schema
from mat3ra.made.material import Material

from ..... import MaterialWithBuildMetadata
from .builder import StrainBuilder
from .configuration import StrainConfiguration


def get_isotropic_strain_matrix(scale_factor: float) -> Matrix3x3Schema:
    """
    Returns a 3x3 isotropic strain matrix for uniform lattice scaling.
    """
    return Matrix3x3Schema(
        root=[
            [scale_factor, 0.0, 0.0],
            [0.0, scale_factor, 0.0],
            [0.0, 0.0, scale_factor],
        ]
    )


def create_strain(
    material: Union[Material, MaterialWithBuildMetadata],
    strain_matrix: Optional[Matrix3x3Schema] = None,
    scale_factor: Optional[float] = None,
) -> MaterialWithBuildMetadata:
    """
    Creates a strained Material based on the provided strain matrix or scale factor.

    Args:
        material (Union[Material, MaterialWithBuildMetadata]): The material to be strained.
        strain_matrix (Optional[Matrix3x3Schema]): The 3x3 strain matrix defining the deformation.
        scale_factor (Optional[float]): An optional scale factor for isotropic scaling instead of a full strain matrix.

    Returns:
        MaterialWithBuildMetadata: The strained material with build metadata.
    """
    if scale_factor is not None:
        if scale_factor <= 0:
            raise ValueError("scale_factor must be positive.")
        strain_matrix = get_isotropic_strain_matrix(scale_factor)
    if strain_matrix is None:
        raise ValueError("Either strain_matrix or scale_factor must be provided.")

    configuration = StrainConfiguration(material=material, strain_matrix=strain_matrix, scale_factor=scale_factor)
    builder = StrainBuilder()
    return builder.get_material(configuration)
