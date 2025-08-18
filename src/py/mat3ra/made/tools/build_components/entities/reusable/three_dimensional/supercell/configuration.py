import numpy as np
from mat3ra.esse.models.materials_category_components.entities.auxiliary.three_dimensional.supercell_matrix_3d import (
    SupercellMatrix3DSchema,
)

from ...base_builder import BaseConfigurationPydantic


class SupercellConfiguration(BaseConfigurationPydantic):
    type: str = "SupercellConfiguration"
    supercell_matrix: SupercellMatrix3DSchema = SupercellMatrix3DSchema(root=np.eye(3).tolist())
