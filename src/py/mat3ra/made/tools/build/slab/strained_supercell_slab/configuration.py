from typing import List, Optional

from mat3ra.esse.models.core.abstract.matrix_3x3 import Matrix3x3Schema
from mat3ra.esse.models.materials_category.pristine_structures.two_dimensional.slab_strained_supercell import (
    SlabStrainedSupercellConfigurationSchema,
)

from mat3ra.made.tools.build.slab.slab.configuration import SlabConfiguration


class SlabStrainedSupercellConfiguration(SlabConfiguration, SlabStrainedSupercellConfigurationSchema):
    strain_matrix: Matrix3x3Schema = SlabStrainedSupercellConfigurationSchema.model_fields["strain_matrix"].default
    xy_supercell_matrix: Optional[List[List[int]]] = SlabStrainedSupercellConfigurationSchema.model_fields[
        "xy_supercell_matrix"
    ].default
