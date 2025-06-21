from typing import Any, List

from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum

from mat3ra.made.material import Material
from mat3ra.made.tools.build import BaseBuilder
from mat3ra.made.tools.build.slab.configuration import (
    AtomicLayersUniqueRepeatedConfiguration,
)
from mat3ra.made.tools.build.stack.configuration import StackConfiguration
from mat3ra.made.tools.build.vacuum.builders import VacuumBuilder
from mat3ra.made.tools.build.vacuum.configuration import VacuumConfiguration
from mat3ra.made.tools.operations.core.binary import stack


class StackBuilder2Components(BaseBuilder):
    _ConfigurationType = StackConfiguration

    def _configuration_to_material(self, configuration_or_material: Any) -> Material:
        if isinstance(configuration_or_material, Material):
            return configuration_or_material
        if isinstance(configuration_or_material, AtomicLayersUniqueRepeatedConfiguration):
            # Local import to avoid circular dependency
            from mat3ra.made.tools.build.slab.builders import AtomicLayersUniqueRepeatedBuilder

            builder = AtomicLayersUniqueRepeatedBuilder()
            return builder.get_material(configuration_or_material)
        if isinstance(configuration_or_material, VacuumConfiguration):
            builder = VacuumBuilder()
            return builder.get_material(configuration_or_material)
        raise ValueError(f"Unknown configuration type: {type(configuration_or_material)}")

    def _generate(self, configuration: StackConfiguration) -> List[Material]:
        first_entity_config = configuration.stack_components[0]
        first_material = self._configuration_to_material(first_entity_config)
        second_entity_config = configuration.stack_components[1]
        second_material = self._configuration_to_material(second_entity_config)

        # Stack the two materials
        stacked_materials = stack([first_material, second_material], configuration.direction or AxisEnum.z)
        return [stacked_materials]
