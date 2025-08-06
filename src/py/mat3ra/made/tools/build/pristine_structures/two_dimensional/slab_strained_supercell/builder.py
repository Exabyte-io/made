from typing import Type

from .configuration import SlabStrainedSupercellConfiguration
from ..slab import SlabBuilder
from .....build_components import MaterialWithBuildMetadata
from .....operations.core.unary import strain, supercell


class SlabStrainedSupercellBuilder(SlabBuilder):
    _ConfigurationType: Type[SlabStrainedSupercellConfiguration] = SlabStrainedSupercellConfiguration

    def _generate(self, configuration: SlabStrainedSupercellConfiguration) -> MaterialWithBuildMetadata:
        slab_material = super()._generate(configuration)
        if configuration.xy_supercell_matrix:
            slab_material = supercell(slab_material, configuration.xy_supercell_matrix)

        strained_slab_material = strain(
            slab_material,
            configuration.strain_matrix,
        )

        return strained_slab_material
