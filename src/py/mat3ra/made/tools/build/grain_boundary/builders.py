from typing import List, Any

from mat3ra.made.material import Material
from .configuration import GrainBoundaryLinearConfiguration, SlabGrainBoundaryConfiguration
from ..interface import ZSLStrainMatchingInterfaceBuilderParameters, InterfaceConfiguration
from ..interface.builders import (
    ZSLStrainMatchingInterfaceBuilder,
    CommensurateLatticeInterfaceBuilderParameters,
)
from ..slab.builders import SlabBuilder, SlabBuilderParameters
from ..slab.builders import SlabWithGapBuilder
from ..slab.configurations import SlabConfiguration
from ..slab.configurations.strained_configurations import SlabStrainedSupercellWithGapConfiguration
from ..stack.builders import Stack2ComponentsBuilder
from ..supercell import create_supercell
from ...analyze.other import get_chemical_formula
from ...third_party import PymatgenInterface


class SlabGrainBoundaryBuilderParameters(ZSLStrainMatchingInterfaceBuilderParameters):
    default_index: int = 0


class GrainBoundaryLinearBuilderParameters(CommensurateLatticeInterfaceBuilderParameters):
    edge_inclusion_tolerance: float = 1.0
    distance_tolerance: float = 1.0


class SlabGrainBoundaryBuilder(ZSLStrainMatchingInterfaceBuilder):
    """
    A builder for creating grain boundaries.

    The grain boundary is created by:
    1. creating an interface between two phases,
    2. then rotating the interface by 90 degrees.
    3. Finally, creating a slab from the rotated interface.
    """

    _BuildParametersType: type(SlabGrainBoundaryBuilderParameters) = SlabGrainBoundaryBuilderParameters  # type: ignore
    _DefaultBuildParameters = SlabGrainBoundaryBuilderParameters()
    _ConfigurationType: type(SlabGrainBoundaryConfiguration) = SlabGrainBoundaryConfiguration  # type: ignore
    _GeneratedItemType: PymatgenInterface = PymatgenInterface  # type: ignore
    selector_parameters: type(  # type: ignore
        SlabGrainBoundaryBuilderParameters
    ) = SlabGrainBoundaryBuilderParameters()  # type: ignore

    def _generate(self, configuration: _ConfigurationType) -> List[_GeneratedItemType]:  # type: ignore
        interface_config = InterfaceConfiguration(
            film_configuration=configuration.phase_1_configuration,
            substrate_configuration=configuration.phase_2_configuration,
            film_termination=configuration.phase_1_termination,
            substrate_termination=configuration.phase_2_termination,
            distance_z=configuration.gap,
            vacuum=configuration.gap,
        )
        return super()._generate(interface_config)

    def _finalize(self, materials: List[Material], configuration: _ConfigurationType) -> List[Material]:
        rot_90_degree_matrix = [[0, 0, 1], [0, 1, 0], [-1, 0, 0]]
        rotated_interfaces = [
            create_supercell(material, supercell_matrix=rot_90_degree_matrix) for material in materials
        ]
        final_slabs: List[Material] = []
        for interface in rotated_interfaces:
            final_slab_config = SlabConfiguration.from_parameters(
                material_or_dict=interface,
                miller_indices=configuration.slab_configuration.miller_indices,
                number_of_layers=configuration.slab_configuration.number_of_layers,
                vacuum=configuration.slab_configuration.vacuum,
            )
            slab_builder_parameters = SlabBuilderParameters(
                xy_supercell_matrix=configuration.slab_configuration.xy_supercell_matrix,
                use_orthogonal_c=True,
            )
            builder = SlabBuilder(build_parameters=slab_builder_parameters)
            final_slab = builder.get_material(final_slab_config)
            final_slabs.append(final_slab)

        return super()._finalize(final_slabs, configuration)

    def _update_material_name(self, material: Material, configuration: _ConfigurationType) -> Material:
        phase_1_formula = get_chemical_formula(configuration.phase_1_configuration.bulk)
        phase_2_formula = get_chemical_formula(configuration.phase_2_configuration.bulk)
        phase_1_miller_indices = "".join([str(i) for i in configuration.phase_1_configuration.miller_indices])
        phase_2_miller_indices = "".join([str(i) for i in configuration.phase_2_configuration.miller_indices])
        new_name = (
            f"{phase_1_formula}({phase_1_miller_indices})-{phase_2_formula}({phase_2_miller_indices}), Grain Boundary"
        )
        material.name = new_name
        return material


class GrainBoundaryLinearBuilder(Stack2ComponentsBuilder):
    """
    Creates a linear grain boundary by stacking two materials along x or y direction.
    """

    _ConfigurationType = GrainBoundaryLinearConfiguration
    _BuildParametersType = GrainBoundaryLinearBuilderParameters
    _DefaultBuildParameters = GrainBoundaryLinearBuilderParameters()

    def _configuration_to_material(self, configuration_or_material: Any) -> Material:
        if isinstance(configuration_or_material, SlabStrainedSupercellWithGapConfiguration):
            builder = SlabWithGapBuilder()
            return builder.get_material(configuration_or_material)
        return super()._configuration_to_material(configuration_or_material)

    def _update_material_name(self, material: Material, configuration: GrainBoundaryLinearConfiguration) -> Material:
        new_name = f"Grain Boundary Linear ({configuration.direction.value})"
        material.name = new_name
        return material
