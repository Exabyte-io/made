from typing import List, Any
from typing import Type

from mat3ra.made.material import Material
from .configuration import GrainBoundaryLinearConfiguration
from ..interface.builders import (
    CommensurateLatticeInterfaceBuilderParameters,
)
from ..slab.builders import SlabBuilder, SlabBuilderParameters
from ..slab.builders import SlabWithGapBuilder
from ..slab.configurations import SlabConfiguration
from ..slab.configurations.strained_configurations import SlabStrainedSupercellWithGapConfiguration
from ..stack.builders import Stack2ComponentsBuilder
from ..supercell import create_supercell
from .configuration import GrainBoundaryConfiguration
from ..interface.builders import InterfaceBuilder, InterfaceBuilderParameters
from ...analyze.other import get_chemical_formula
from ...modify import wrap_to_unit_cell
from ...operations.core.unary import supercell, edit_cell
from ...utils import AXIS_TO_INDEX_MAP


class GrainBoundaryBuilderParameters(InterfaceBuilderParameters):
    pass


class GrainBoundaryLinearBuilderParameters(CommensurateLatticeInterfaceBuilderParameters):
    edge_inclusion_tolerance: float = 1.0
    distance_tolerance: float = 1.0


class GrainBoundaryBuilder(InterfaceBuilder):
    """
    Builder for creating grain boundaries.

    Uses InterfaceBuilder to create a Z-stacked interface, then transforms the supercell
    to flip the Z direction to X direction for grain boundary orientation.
    """

    _BuildParametersType: Type[GrainBoundaryBuilderParameters] = GrainBoundaryBuilderParameters
    _DefaultBuildParameters = GrainBoundaryBuilderParameters()

    def _generate(self, configuration: GrainBoundaryConfiguration) -> Material:
        interface = super()._generate(configuration)
        rotated_interface = supercell(interface, [[0, 0, 1], [0, 1, 0], [1, 0, 0]])
        wrapped_interface = wrap_to_unit_cell(rotated_interface)
        return wrapped_interface

    def _update_material_name(self, material: Material, configuration: GrainBoundaryConfiguration) -> Material:
        phase_1_formula = get_chemical_formula(configuration.phase_1_configuration.atomic_layers.crystal)
        phase_2_formula = get_chemical_formula(configuration.phase_2_configuration.atomic_layers.crystal)
        phase_1_miller = "".join([str(i) for i in configuration.phase_1_configuration.atomic_layers.miller_indices])
        phase_2_miller = "".join([str(i) for i in configuration.phase_2_configuration.atomic_layers.miller_indices])
        new_name = f"{phase_1_formula}({phase_1_miller})-{phase_2_formula}({phase_2_miller}), Grain Boundary"
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

    def _generate(self, configuration: GrainBoundaryLinearConfiguration) -> Material:
        stacked_material = super()._generate(configuration)
        return stacked_material

    def _get_axis_index(self, direction_str):
        return AXIS_TO_INDEX_MAP[direction_str]
