from typing import Type

from mat3ra.made.tools.build import MaterialWithBuildMetadata
from mat3ra.made.tools.build.slab.atomic_layers_unique_repeated.builder import AtomicLayersUniqueRepeatedBuilder
from mat3ra.made.tools.build.slab.slab.builder_parameters import SlabBuilderParameters
from mat3ra.made.tools.build.slab.slab.configuration import SlabConfiguration
from mat3ra.made.tools.build.slab.utils import get_orthogonal_c_slab
from mat3ra.made.tools.build.stack.builder import StackNComponentsBuilder
from mat3ra.made.tools.operations.core.unary import supercell


class SlabBuilder(StackNComponentsBuilder):
    _BuildParametersType: Type[SlabBuilderParameters] = SlabBuilderParameters
    _DefaultBuildParameters: SlabBuilderParameters = SlabBuilderParameters()

    @property
    def stack_component_types_conversion_map(self):
        return {
            **super().stack_component_types_conversion_map,
            SlabConfiguration: AtomicLayersUniqueRepeatedBuilder,
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
