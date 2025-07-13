import numpy as np
from mat3ra.esse.models.core.abstract.matrix_3x3 import Matrix3x3Schema
from mat3ra.esse.models.materials_category.pristine_structures.two_dimensional.slab_strained_supercell import (
    SlabStrainedSupercellConfigurationSchema,
)

from .slab_configuration import SlabConfiguration


class SlabStrainedSupercellConfiguration(SlabConfiguration, SlabStrainedSupercellConfigurationSchema):
    type: str = "SlabStrainedSupercellConfiguration"
    strain_matrix: Matrix3x3Schema = Matrix3x3Schema(root=np.eye(3).tolist())
