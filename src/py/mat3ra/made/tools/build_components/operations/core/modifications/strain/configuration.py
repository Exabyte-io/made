from typing import Union

from mat3ra.esse.models.core.abstract.matrix_3x3 import Matrix3x3Schema
from mat3ra.made.material import Material
from mat3ra.made.tools.build_components.entities.reusable.base_builder import BaseConfigurationPydantic

from ..... import MaterialWithBuildMetadata


class StrainConfiguration(BaseConfigurationPydantic):
    type: str = "StrainConfiguration"
    material: Union[Material, MaterialWithBuildMetadata]
    strain_matrix: Matrix3x3Schema
