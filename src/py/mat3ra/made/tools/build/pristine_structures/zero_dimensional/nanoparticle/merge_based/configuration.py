from typing import List

from mat3ra.esse.models.materials_category_components.operations.core.combinations.merge import MergeMethodsEnum

from mat3ra.made.tools.build_components.entities.auxiliary.zero_dimensional.void_region.configuration import \
    VoidRegionConfiguration
from mat3ra.made.tools.build_components.entities.reusable.two_dimensional.atomic_layers.configuration import \
    SlabConfiguration
from mat3ra.made.tools.build_components.operations.core.combinations.merge.configuration import MergeConfiguration


class NanoparticleConfiguration(MergeConfiguration):
    merge_components: List = [SlabConfiguration, VoidRegionConfiguration]
    merge_method: MergeMethodsEnum = MergeMethodsEnum.REPLACE

    @property
    def void_region_configuration(self) -> VoidRegionConfiguration:
        return self.merge_components[1]
