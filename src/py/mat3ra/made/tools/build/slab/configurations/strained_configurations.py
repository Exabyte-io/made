from typing import Optional

import numpy as np
from mat3ra.esse.models.core.abstract.matrix_3x3 import Matrix3x3Schema
from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from mat3ra.esse.models.materials_category.pristine_structures.two_dimensional.slab_strained_supercell import (
    SlabStrainedSupercellConfigurationSchema,
)
from mat3ra.esse.models.materials_category.pristine_structures.two_dimensional.slab_strained_supercell_with_gap import (
    SlabStrainedSupercellWithGapConfigurationSchema,
)

from .slab_configuration import SlabConfiguration
from ...supercell import SupercellConfiguration


class SlabStrainedSupercellConfiguration(SlabConfiguration, SlabStrainedSupercellConfigurationSchema):
    type: str = "SlabStrainedSupercellConfiguration"
    strain_matrix: Matrix3x3Schema = Matrix3x3Schema(root=np.eye(3).tolist())


class SlabStrainedSupercellWithGapConfiguration(
    SlabStrainedSupercellConfiguration, SlabStrainedSupercellWithGapConfigurationSchema
):
    type: str = "SlabStrainedSupercellWithGapConfiguration"
    gap: Optional[float] = None  # If provided, the film is shifted to have it as smallest distance to the top.
    gap_direction: AxisEnum = AxisEnum.z
