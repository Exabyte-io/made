from typing import List

from mat3ra.esse.models.materials_category_components.operations.core.combinations.merge import MergeMethodsEnum

from ...merge.configuration import MergeConfiguration
from ...slab.slab.configuration import SlabConfiguration
from ...void_region.configuration import VoidRegionConfiguration


class NanoparticleConfiguration(MergeConfiguration):
    merge_components: List = [SlabConfiguration, VoidRegionConfiguration]
    merge_method: MergeMethodsEnum = MergeMethodsEnum.REPLACE

    @property
    def void_region_configuration(self) -> VoidRegionConfiguration:
        return self.merge_components[1]