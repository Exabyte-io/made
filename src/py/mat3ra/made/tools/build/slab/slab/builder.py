from typing import Type

from ... import MaterialWithBuildMetadata
from ..atomic_layers_unique_repeated.builder import AtomicLayersUniqueRepeatedBuilder
from ..atomic_layers_unique_repeated.configuration import (
    AtomicLayersUniqueRepeatedConfiguration,
)
from .build_parameters import SlabBuilderParameters
from .configuration import SlabConfiguration
from ..utils import get_orthogonal_c_slab
from ...stack.builder import StackNComponentsBuilder
from ....operations.core.unary import supercell


class SlabBuilder(StackNComponentsBuilder):
    _BuildParametersType: Type[SlabBuilderParameters] = SlabBuilderParameters
    _DefaultBuildParameters: SlabBuilderParameters = SlabBuilderParameters()

    @property
    def stack_component_types_conversion_map(self):
        return {
            **super().stack_component_types_conversion_map,
            AtomicLayersUniqueRepeatedConfiguration: AtomicLayersUniqueRepeatedBuilder,
        }

    def _generate(self, configuration: SlabConfiguration) -> MaterialWithBuildMetadata:
        stack_as_material = super()._generate(configuration)
        supercell_slab = supercell(stack_as_material, self.build_parameters.xy_supercell_matrix)
        if self.build_parameters.use_orthogonal_c:
            supercell_slab = get_orthogonal_c_slab(supercell_slab)
        return supercell_slab

    def _update_material_name(
        self, material: MaterialWithBuildMetadata, configuration: SlabConfiguration
    ) -> MaterialWithBuildMetadata:
        # for example: "Si(001), termination Si_P4/mmm_1, Slab"
        material = AtomicLayersUniqueRepeatedBuilder()._update_material_name(material, configuration.atomic_layers)
        new_name = f"{material.name}, Slab"
        material.name = new_name
        return material
