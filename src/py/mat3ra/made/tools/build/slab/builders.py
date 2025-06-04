from typing import Any

from mat3ra.esse.models.apse.materials.builders.slab.pymatgen.parameters import (
    PymatgenSlabGeneratorParametersSchema,
)
from pydantic import BaseModel

from mat3ra.made.material import Material
from .configuration import (
    SlabConfiguration,
    AtomicLayersUniqueRepeatedConfiguration,
    VacuumConfiguration,
    AtomicLayersUniqueRepeatedBuilder,
    StackConfiguration,
)
from .termination import Termination
from ..vacuum.builders import VacuumBuilder
from ...build import BaseBuilder
from ...operations.core.unary import stack, supercell


class SlabSelectorParameters(BaseModel):
    termination: Termination


class PymatgenSlabGeneratorParameters(PymatgenSlabGeneratorParametersSchema):
    # Parameters described in https://github.com/materialsproject/pymatgen/blob/585bb673c4aa222669c4b0d72ffeec3dbf092630/pymatgen/core/surface.py#L1187
    pass


class SlabBuilderParameters(PymatgenSlabGeneratorParameters):
    make_primitive: bool = False
    use_orthogonal_c: bool = True


class StackBuilder2Components(BaseBuilder):

    def configuration_to_material(self, configuration_or_material: Any) -> Material:
        if isinstance(configuration_or_material, Material):
            return configuration_or_material
        if isinstance(configuration_or_material, AtomicLayersUniqueRepeatedConfiguration):
            builder = AtomicLayersUniqueRepeatedBuilder()
            return builder.get_material(configuration_or_material)
        if isinstance(configuration_or_material, VacuumConfiguration):
            builder = VacuumBuilder()
            return builder.get_material(configuration_or_material)
        # If we reach here, we don't know how to handle this configuration
        raise ValueError(f"Unknown configuration type: {type(configuration_or_material)}")

    def generate(self, configuration: StackConfiguration) -> Material:
        first_entity_config = configuration.stack_components[0]
        first_material = self.configuration_to_material(first_entity_config)
        second_entity_config = configuration.stack_components[1]
        second_material = self.configuration_to_material(second_entity_config)

        # Stack the two materials
        stacked_materials = stack([first_material, second_material], "z")
        return stacked_materials

    def get_material(self, configuration: SlabConfiguration) -> Material:
        return self.generate(configuration)


class SlabBuilder(StackBuilder2Components):

    def generate(self, configuration: SlabConfiguration) -> Material:
        # Use builders to construct materials from configurations
        atomic_layers_config = configuration.atomic_layers
        vacuum_config = configuration.vacuum_configuration

        # Build the atomic layers material
        if isinstance(atomic_layers_config, AtomicLayersUniqueRepeatedConfiguration):
            atomic_builder = AtomicLayersUniqueRepeatedBuilder()
            repeated_layers = atomic_builder.get_material(atomic_layers_config)
        else:
            # For AtomicLayersUnique, just use the orthogonal_c_cell
            repeated_layers = atomic_layers_config.orthogonal_c_cell

        # Build the vacuum material
        vacuum_builder = VacuumBuilder()
        vacuum_layer = vacuum_builder.get_material(vacuum_config)

        # Stack the materials
        stacked_materials = stack([repeated_layers, vacuum_layer], "z")
        supercell_slab = supercell(stacked_materials, configuration.xy_supercell_matrix)
        return supercell_slab

    def get_material(self, configuration: SlabConfiguration) -> Material:
        return self.generate(configuration)


#
# class SlabBuilder(ConvertGeneratedItemsPymatgenStructureMixin, BaseBuilder):
#     build_parameters: Optional[SlabBuilderParameters]
#     _DefaultBuildParameters: SlabBuilderParameters = SlabBuilderParameters()
#     _ConfigurationType: type(SlabConfiguration) = SlabConfiguration  # type: ignore
#     _GeneratedItemType: Material = Material  # type: ignore
#     _SelectorParametersType: type(SlabSelectorParameters) = SlabSelectorParameters  # type: ignore
#
#     def get_material(self, configuration: _ConfigurationType) -> Material:  # type: ignore
#         self._configuration = configuration
#         return super().get_material(configuration)
#
#     def _generate(self, configuration: _ConfigurationType) -> List[Material]:  # type: ignore
#         atomic_layers: AtomicLayersUniqueRepeatedConfiguration = configuration.atomic_layers
#         vacuum: VacuumConfiguration = configuration.stack_components[1]
#         params = self.build_parameters or self._DefaultBuildParameters
#
#         slab_materials = atomic_layers.get_slabs(
#             min_slab_size=self._configuration.number_of_layers,
#             min_vacuum_size=params.min_vacuum_size,
#             in_unit_planes=params.in_unit_planes,
#             make_primitive=params.make_primitive,
#             symmetrize=params.symmetrize,
#             use_orthogonal_c=params.use_orthogonal_c,
#         )
#
#         stacked_materials = []
#         for slab_material in slab_materials:
#             stacked = stack_two_components(slab_material, vacuum, direction=configuration.direction)
#             stacked_materials.append(stacked)
#
#         supercells_of_stacked_materials = [
#             create_supercell(material, configuration.xy_supercell_matrix) for material in stacked_materials
#         ]
#
#         return supercells_of_stacked_materials
#
#     def _update_material_name(self, material: Material, configuration: SlabConfiguration) -> Material:
#         atomic_layers = configuration.atomic_layers
#
#         formula = get_chemical_formula(atomic_layers.crystal)
#         miller_indices = "".join(str(i) for i in atomic_layers.miller_indices)
#         termination = configuration.termination
#         material.name = f"{formula}({miller_indices}), termination {termination}, Slab"
#         return material
