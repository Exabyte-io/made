from typing import Optional, Union

from mat3ra.esse.models.core.abstract.matrix_3x3 import Matrix3x3Schema
from mat3ra.made.material import Material
from mat3ra.made.tools.build_components.entities.reusable.base_builder import BaseConfigurationPydantic

from ..... import MaterialWithBuildMetadata


class StrainConfiguration(BaseConfigurationPydantic):
    """
    Configuration for a strain operation on a Material.

    Args:
        material (Material): The Material object to be strained.
        strain_matrix (Matrix3x3Schema): The 3x3 strain matrix defining the deformation.
        scale_factor (Optional[float]): An optional scale factor for isotropic scaling instead of a full strain matrix.
    """

    type: str = "StrainConfiguration"
    material: Union[Material, MaterialWithBuildMetadata]
    strain_matrix: Matrix3x3Schema
    scale_factor: Optional[float] = None
