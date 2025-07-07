from typing import List, Union

from mat3ra.esse.models.materials_category_components.operations.core.combinations.merge import MergeMethodsEnum
from mat3ra.made.material import Material
from mat3ra.made.tools.build.merge.configuration import MergeConfiguration
from mat3ra.made.tools.build import BaseConfigurationPydantic


class SlabDefectConfiguration(MergeConfiguration):
    """
    Configuration for merging a slab with additional layers and an isolated defect.

    Args:
        merge_components: List containing [slab, isolated_defect].
        merge_method: Method to use for merging (default: add).
        auto_add_vacuum: Whether to automatically add vacuum if needed.
        vacuum_thickness: Thickness of vacuum to add if auto_add_vacuum is True.
    """

    type: str = "SlabDefectConfiguration"
    auto_add_vacuum: bool = True
    vacuum_thickness: float = 5.0
    merge_method = MergeMethodsEnum.ADD

    @property
    def slab(self) -> Material:
        return self.merge_components[0]

    @property
    def isolated_defect(self) -> Material:
        return self.merge_components[1]

    @classmethod
    def from_materials(
        cls,
        slab: Material,
        isolated_defect: Material,
        merge_method: MergeMethodsEnum = MergeMethodsEnum.add,
        auto_add_vacuum: bool = True,
        vacuum_thickness: float = 5.0,
    ) -> "SlabDefectConfiguration":
        """
        Creates a SlabDefectConfiguration from materials.

        Args:
            slab: The slab material with additional layers.
            isolated_defect: The isolated defect material.
            merge_method: Method to use for merging.
            auto_add_vacuum: Whether to automatically add vacuum if needed.
            vacuum_thickness: Thickness of vacuum to add if auto_add_vacuum is True.

        Returns:
            SlabDefectConfiguration: The merge configuration.
        """
        return cls(
            merge_components=[slab, isolated_defect],
            merge_method=merge_method,
            auto_add_vacuum=auto_add_vacuum,
            vacuum_thickness=vacuum_thickness,
        )
