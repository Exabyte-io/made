from typing import Type

from mat3ra.made.tools.build import TConfiguration, MaterialWithBuildMetadata
from mat3ra.made.tools.build.slab.slab.builder import SlabBuilder
from mat3ra.made.tools.build.slab.strained_supercell_slab.configuration import SlabStrainedSupercellConfiguration
from mat3ra.made.tools.operations.core.unary import supercell, strain


class SlabStrainedSupercellBuilder(SlabBuilder):
    _ConfigurationType: Type[SlabStrainedSupercellConfiguration] = SlabStrainedSupercellConfiguration

    def _generate(self, configuration: TConfiguration) -> MaterialWithBuildMetadata:
        slab_material = super()._generate(configuration)
        if configuration.xy_supercell_matrix:
            slab_material = supercell(slab_material, configuration.xy_supercell_matrix)

        strained_slab_material = strain(
            slab_material,
            configuration.strain_matrix,
        )

        return strained_slab_material
