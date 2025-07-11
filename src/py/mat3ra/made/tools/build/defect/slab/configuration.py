from mat3ra.made.material import Material
from mat3ra.made.tools.build.stack.configuration import StackConfiguration
from mat3ra.made.tools.build.vacuum.configuration import VacuumConfiguration


class SlabStackConfiguration(StackConfiguration):
    """
    Configuration for stacking a slab, an isolated defect and a vacuum layer.

    Args:
        stack_components: List containing [slab, slab_component, vacuum].
        direction: Direction of stacking (default: z).
    """

    type: str = "SlabStackConfiguration"

    @property
    def slab(self) -> Material:
        return self.stack_components[0]

    @property
    def added_component(self) -> Material:
        return self.stack_components[1]

    @property
    def vacuum_configuration(self) -> VacuumConfiguration:
        return self.stack_components[2]


class IslandDefectConfiguration(SlabStackConfiguration):
    """
    Configuration for creating an island defect on a slab surface.

    Args:
        stack_components: List containing [slab, isolated_island, vacuum].
    """

    type: str = "IslandDefectConfiguration"
