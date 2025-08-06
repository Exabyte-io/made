from typing import Type

from .configuration import SlabStrainedSupercellConfiguration
from .....analyze.build_metadata_analyzer import TypeConfiguration
from .....build_components import MaterialWithBuildMetadata
from .....build_components.entities.reusable.two_dimensional.atomic_layers.builder import SlabBuilder
from .....operations.core.unary import strain, supercell


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
