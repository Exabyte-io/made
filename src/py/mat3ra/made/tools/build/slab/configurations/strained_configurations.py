from typing import Optional

from mat3ra.esse.models.core.abstract.matrix_3x3 import Matrix3x3Schema
from mat3ra.esse.models.materials_category.pristine_structures.two_dimensional.slab_strained_supercell import (
    SlabStrainedSupercellConfigurationSchema,
)
from mat3ra.esse.models.materials_category.pristine_structures.two_dimensional.slab_strained_supercell_with_gap import (
    SlabStrainedSupercellWithGapConfigurationSchema,
)

from .slab_configuration import SlabConfiguration


class SlabStrainedSupercellConfiguration(SlabConfiguration, SlabStrainedSupercellConfigurationSchema):
    type: str = "SlabStrainedSupercellConfiguration"
    strain_matrix: Matrix3x3Schema = SlabStrainedSupercellConfigurationSchema.model_fields["strain_matrix"].default


class SlabStrainedSupercellWithGapConfiguration(
    SlabStrainedSupercellConfiguration, SlabStrainedSupercellWithGapConfigurationSchema
):
    type: str = "SlabStrainedSupercellWithGapConfiguration"
    gap: Optional[float] = None  # If provided, the film is shifted to have it as smallest distance to the top.
