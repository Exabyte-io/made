from typing import Type

from mat3ra.made.tools.build_components import MaterialWithBuildMetadata
from mat3ra.made.tools.operations.core.unary import strain, supercell

from ..slab.builder import SlabBuilder
from .configuration import SlabStrainedSupercellConfiguration



class SlabStrainedSupercellBuilder(SlabBuilder):
    _ConfigurationType: Type[SlabStrainedSupercellConfiguration] = SlabStrainedSupercellConfiguration

    def _generate(self, configuration: TypeConfiguration) -> MaterialWithBuildMetadata:
        slab_material = super()._generate(configuration)
        if configuration.xy_supercell_matrix:
            slab_material = supercell(slab_material, configuration.xy_supercell_matrix)

        strained_slab_material = strain(
            slab_material,
            configuration.strain_matrix,
        )

        return strained_slab_material
