from typing import Dict, Optional, List

from mat3ra.esse.models.materials_category_components.entities.auxiliary.zero_dimensional.void_region import (
    VoidRegionSchema,
)
from mat3ra.esse.models.materials_category_components.operations.core.combinations.merge import MergeMethodsEnum

from .enums import ASENanoparticleShapesEnum
from ..merge.configuration import MergeConfiguration
from ..slab.configurations import SlabConfiguration


class NanoparticleConfiguration(MergeConfiguration):
    merge_components: List = [SlabConfiguration, VoidRegionSchema]
    merge_method: MergeMethodsEnum = MergeMethodsEnum.ADD


class ASEBasedNanoparticleConfiguration(NanoparticleConfiguration):
    """
    Configuration for building a nanoparticle.

    Attributes:
        shape (NanoparticleShapes): The desired shape of the nanoparticle.
        parameters (dict): Dictionary of parameters to pass to the corresponding ASE constructor.
    """

    shape: ASENanoparticleShapesEnum
    parameters: Optional[Dict] = None  # Shape-specific parameters (e.g., layers, size)

    @property
    def element(self) -> str:
        return self.material.basis.elements.get_element_value_by_index(0)

    @property
    def _json(self):
        return {
            **super()._json,
            "shape": self.shape.value,
            "parameters": self.parameters,
        }
