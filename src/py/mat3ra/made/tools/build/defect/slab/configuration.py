from mat3ra.esse.models.materials_category_components.operations.core.combinations.merge import MergeMethodsEnum

from mat3ra.made.material import Material
from mat3ra.made.tools.build.merge.configuration import MergeConfiguration
from mat3ra.made.tools.build.stack.configuration import StackConfiguration
from mat3ra.made.tools.build.vacuum.configuration import VacuumConfiguration


# fmt: off
# fmt: on


class SlabDefectConfiguration(MergeConfiguration):
    """
    Configuration for merging a slab with additional layers and an isolated defect.

    Args:
        merge_components: List containing [slab, isolated_defect].
        merge_method: Method to use for merging (default: add).
    """

    type: str = "SlabDefectConfiguration"
    merge_method: MergeMethodsEnum = MergeMethodsEnum.ADD

    @property
    def slab(self) -> Material:
        return self.merge_components[0]

    @property
    def isolated_defect(self) -> Material:
        return self.merge_components[1]


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
