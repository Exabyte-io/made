from typing import Union

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
    strain_matrix: Matrix3x3Schema,
) -> MaterialWithBuildMetadata:
    configuration = StrainConfiguration(material=material, strain_matrix=strain_matrix)
    builder = StrainBuilder()
    return builder.get_material(configuration)
