from typing import Optional

import numpy as np
from mat3ra.esse.models.core.abstract.matrix_3x3 import Matrix3x3Schema
from mat3ra.esse.models.materials_category.pristine_structures.two_dimensional.slab_strained_supercell import (
    SlabStrainedSupercellConfigurationSchema,
)
from mat3ra.esse.models.materials_category_components.entities.auxiliary.two_dimensional.supercell_matrix_2d import (
    SupercellMatrix2DSchema,
)

from ..slab.configuration import SlabConfiguration


class SlabStrainedSupercellConfiguration(SlabConfiguration, SlabStrainedSupercellConfigurationSchema):
    strain_matrix: Matrix3x3Schema = Matrix3x3Schema(root=np.eye(3).tolist())
    xy_supercell_matrix: Optional[SupercellMatrix2DSchema] = SupercellMatrix2DSchema(root=[[1, 0], [0, 1]])
