from typing import Any, TypeVar, Optional

from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum

from mat3ra.made.material import Material
from mat3ra.made.tools.build import BaseSingleBuilder
from mat3ra.made.tools.build.stack.configuration import StackConfiguration
from mat3ra.made.tools.build.vacuum.builders import VacuumBuilder
from mat3ra.made.tools.build.vacuum.configuration import VacuumConfiguration
from mat3ra.made.tools.operations.core.binary import stack


class Stack2ComponentsBuilder(BaseSingleBuilder):
    _ConfigurationType = StackConfiguration

    def _set_vacuum_crystal(self, stack_components):
        previous_material = None
        for entity_config in stack_components:
            if isinstance(entity_config, VacuumConfiguration) and entity_config.crystal is None:
                if previous_material is None:
                    raise ValueError(
                        "VacuumConfiguration.crystal is None and no previous material to use as reference."
                    )
                entity_config.crystal = previous_material
            if isinstance(entity_config, Material):
                previous_material = entity_config
            else:
                previous_material = self._configuration_to_material(entity_config)

    def _configuration_to_material(self, configuration_or_material: Any) -> Optional[Material]:
        if configuration_or_material is None:
            return None  # Return None for None input - caller should handle appropriately
        if isinstance(configuration_or_material, Material):
            return configuration_or_material
        if isinstance(configuration_or_material, VacuumConfiguration):
            builder = VacuumBuilder()
            return builder.get_material(configuration_or_material)
        raise ValueError(f"Unknown configuration type: {type(configuration_or_material)}")

    def _generate(self, configuration: StackConfiguration) -> Material:
        self._set_vacuum_crystal(configuration.stack_components)
        first_entity_config = configuration.stack_components[0]
        first_material = self._configuration_to_material(first_entity_config)
        second_entity_config = configuration.stack_components[1]
        second_material = self._configuration_to_material(second_entity_config)
        stacked_material = stack([first_material, second_material], configuration.direction or AxisEnum.z)
        return stacked_material


StackConfigurationType = TypeVar("StackConfigurationType", bound=StackConfiguration)


class StackNComponentsBuilder(Stack2ComponentsBuilder):
    _ConfigurationType = StackConfiguration

    def _generate(self, configuration: StackConfigurationType) -> Material:
        self._set_vacuum_crystal(configuration.stack_components)
        materials = []
        for entity_config in configuration.stack_components:
            material = self._configuration_to_material(entity_config)
            if material is not None:
                materials.append(material)
        stacked_material = stack(materials, configuration.direction or AxisEnum.z)
        return stacked_material
