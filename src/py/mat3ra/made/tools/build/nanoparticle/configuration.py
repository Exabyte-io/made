from typing import Dict, Optional, List

from mat3ra.esse.models.materials_category_components.entities.auxiliary.zero_dimensional.void_region import (
    VoidRegionSchema,
)
from mat3ra.esse.models.materials_category_components.operations.core.combinations.merge import MergeMethodsEnum

from .enums import NanoparticleShapesEnum
from .. import BaseConfigurationPydantic
from ..merge.configuration import MergeConfiguration
from ..slab.configurations import SlabConfiguration


class NanoparticleConfiguration(MergeConfiguration):
    merge_components: List = [SlabConfiguration, VoidRegionSchema]
    merge_method: MergeMethodsEnum = MergeMethodsEnum.ADD


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
