from typing import Generic, Tuple, TypeVar

from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from mat3ra.made.tools.analyze.lattice import LatticeMaterialAnalyzer
from mat3ra.made.tools.build import BaseConfigurationPydantic
from mat3ra.made.tools.build.vacuum.configuration import VacuumConfiguration

T = TypeVar("T", bound=BaseConfigurationPydantic)


class BaseStackAnalyzer(LatticeMaterialAnalyzer, Generic[T]):
    """
    Base analyzer for creating stack configurations with a component and vacuum.
    """

    miller_indices_uv: Tuple[int, int]

    def get_component_configuration(self, **kwargs) -> T:
        raise NotImplementedError

    def get_component_builder(self):
        raise NotImplementedError

    def get_stack_direction(self) -> AxisEnum:
        raise NotImplementedError

    def get_vacuum_direction(self) -> AxisEnum:
        raise NotImplementedError

    def get_configuration(self, vacuum_size: float, **kwargs) -> BaseConfigurationPydantic:
        """
        Create a stack configuration with component and vacuum.
        Args:
            vacuum_size: The size of the vacuum region in Angstroms (cartesian).
            **kwargs: Additional arguments passed to get_component_configuration.
        """
        component_config = self.get_component_configuration(**kwargs)

        component_builder = self.get_component_builder()
        component_material = component_builder.get_material(component_config)

        vacuum_config = VacuumConfiguration(
            size=vacuum_size,
            crystal=component_material,
            direction=self.get_vacuum_direction(),
        )

        from mat3ra.made.tools.build.stack.configuration import StackConfiguration

        return StackConfiguration(
            stack_components=[component_config, vacuum_config],
            direction=self.get_stack_direction(),
        )
