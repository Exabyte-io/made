from typing import List
from mat3ra.esse.models.materials_category_components.operations.core.combinations.merge import MergeMethodsEnum
# fmt: off
from mat3ra.esse.models.materials_category.defective_structures.two_dimensional. \
    adatom.configuration import AdatomPointDefectSchema

# fmt: on

from mat3ra.made.material import Material
from mat3ra.made.tools.build.merge.configuration import MergeConfiguration
from mat3ra.made.tools.site import CrystalSite
from mat3ra.made.tools.utils.coordinate import CoordinateCondition


class VoidSite(CrystalSite):
    """
    A void site represents an empty space defined by a coordinate condition.
    This is used to create islands by removing atoms that don't satisfy the condition.
    """
    condition: CoordinateCondition

    def condition_function(self, coordinate: List[float]) -> bool:
        return self.condition.condition(coordinate)


class IslandDefectConfiguration(MergeConfiguration):
    """
    Configuration for creating an island defect on a slab surface.
    
    Args:
        merge_components: List containing [slab, void_site].
        merge_method: Method to use for merging (default: replace).
    """
    
    type: str = "IslandDefectConfiguration"
    merge_method: MergeMethodsEnum = MergeMethodsEnum.REPLACE
    
    @property
    def slab(self) -> Material:
        return self.merge_components[0]
    
    @property
    def void_site(self) -> VoidSite:
        return self.merge_components[1]


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


class AdatomDefectConfiguration(MergeConfiguration, AdatomPointDefectSchema):
    """
    Configuration for creating an adatom defect on a slab surface.S

    Args:
        merge_components: List containing [slab, isolated_defect].
    """

    type: str = "AdatomDefectConfiguration"
