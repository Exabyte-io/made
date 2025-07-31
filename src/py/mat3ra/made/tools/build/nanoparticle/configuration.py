from typing import Dict, Optional, List

from mat3ra.esse.models.materials_category_components.operations.core.combinations.merge import MergeMethodsEnum

from .enums import NanoparticleShapesEnum
from .. import BaseConfigurationPydantic
from ..merge.configuration import MergeConfiguration
from ..slab.configurations import SlabConfiguration
from ..void_region.configuration import VoidRegionConfiguration


class NanoparticleConfiguration(MergeConfiguration):
    merge_components: List = [SlabConfiguration, VoidRegionConfiguration]
    merge_method: MergeMethodsEnum = MergeMethodsEnum.REPLACE

    @property
    def void_region_configuration(self) -> VoidRegionConfiguration:
        return self.merge_components[1]


class ASEBasedNanoparticleConfiguration(BaseConfigurationPydantic):
    """
    Configuration for building a nanoparticle with bottom-up approach from atoms.

    Attributes:
        shape (NanoparticleShapes): The desired shape of the nanoparticle.
        parameters (dict): Dictionary of parameters to pass to the corresponding ASE constructor.
    """

    shape: NanoparticleShapesEnum
    parameters: Optional[Dict] = None
    element: str
    vacuum_padding: float = 10.0
