from typing import Any

from mat3ra.made.material import Material
from mat3ra.made.tools.build.defect.slab.configuration import SlabStackConfiguration
from mat3ra.made.tools.build.slab.builders import SlabBuilder
from mat3ra.made.tools.build.slab.configurations import SlabConfiguration
from mat3ra.made.tools.build.stack.builders import StackNComponentsBuilder


class SlabStackBuilder(StackNComponentsBuilder):
    def _configuration_to_material(self, configuration_or_material: Any) -> Material:
        if isinstance(configuration_or_material, SlabConfiguration):
            return SlabBuilder().get_material(configuration_or_material)
        else:
            return super()._configuration_to_material(configuration_or_material)

    def _generate(self, configuration: SlabStackConfiguration) -> Material:
        return super()._generate(configuration)
