from typing import Any

from mat3ra.made.material import Material
from ..slab.configuration import SlabStackConfiguration
from ...slab.builders import SlabBuilder, SlabWithGapBuilder
from ...slab.configurations import SlabConfiguration, SlabStrainedSupercellWithGapConfiguration
from ...stack.builders import StackNComponentsBuilder


class SlabStackBuilder(StackNComponentsBuilder):
    def _configuration_to_material(self, configuration_or_material: Any) -> Material:
        if isinstance(configuration_or_material, SlabStrainedSupercellWithGapConfiguration):
            return SlabWithGapBuilder().get_material(configuration_or_material)
        if isinstance(configuration_or_material, SlabConfiguration):
            return SlabBuilder().get_material(configuration_or_material)
        else:
            return super()._configuration_to_material(configuration_or_material)

    def _generate(self, configuration: SlabStackConfiguration) -> Material:
        return super()._generate(configuration)
