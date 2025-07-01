from typing import List, Type

from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum

from mat3ra.made.material import Material
from mat3ra.made.tools.build import BaseBuilder, BaseBuilderParameters
from mat3ra.made.tools.build.interface.builders import InterfaceBuilder, InterfaceBuilderParameters
from mat3ra.made.tools.build.slab.builders import SlabStrainedSupercellBuilder
from mat3ra.made.tools.build.stack.builders import Stack2ComponentsBuilder
from mat3ra.made.tools.build.stack.configuration import StackConfiguration
from mat3ra.made.tools.build.utils import stack_two_materials_xy
from mat3ra.made.tools.build.vacuum.configuration import VacuumConfiguration
from mat3ra.made.tools.modify import translate_by_vector
from .configuration import GrainBoundaryConfiguration, GrainBoundaryWithVacuumConfiguration
from ..interface import InterfaceConfiguration
from ...analyze.other import get_chemical_formula


class GrainBoundaryBuilderParameters(InterfaceBuilderParameters):
    pass


class GrainBoundaryBuilder(InterfaceBuilder):
    """
    Builder for creating grain boundaries.

    Inherits from InterfaceBuilder but stacks materials in x direction instead of z,
    with shift applied in y,z directions and gap in x direction.
    """

    _BuildParametersType: Type[GrainBoundaryBuilderParameters] = GrainBoundaryBuilderParameters
    _DefaultBuildParameters = GrainBoundaryBuilderParameters()

    def _generate(self, configuration: GrainBoundaryConfiguration) -> Material:
        interface_config = InterfaceConfiguration(...)

        return super()._generate(interface_config)

    def _update_material_name(self, material: Material, configuration: GrainBoundaryConfiguration) -> Material:
        phase_1_formula = get_chemical_formula(
            configuration.phase_1_configuration.substrate_configuration.atomic_layers.crystal
        )
        phase_2_formula = get_chemical_formula(
            configuration.phase_2_configuration.film_configuration.atomic_layers.crystal
        )
        phase_1_miller = "".join(
            [str(i) for i in configuration.phase_1_configuration.substrate_configuration.atomic_layers.miller_indices]
        )
        phase_2_miller = "".join(
            [str(i) for i in configuration.phase_2_configuration.film_configuration.atomic_layers.miller_indices]
        )

        new_name = f"{phase_1_formula}({phase_1_miller})-{phase_2_formula}({phase_2_miller}), Grain Boundary"
        material.name = new_name
        return material


class GrainBoundaryWithVacuumBuilderParameters(BaseBuilderParameters):
    use_orthogonal_c: bool = True


class GrainBoundaryWithVacuumBuilder(BaseBuilder):
    """
    Builder for creating grain boundaries with vacuum (slab).

    Adds vacuum directly to the grain boundary material instead of trying to create
    a slab from it, since grain boundary materials don't have conventional terminations.
    """

    _BuildParametersType: Type[GrainBoundaryWithVacuumBuilderParameters] = GrainBoundaryWithVacuumBuilderParameters
    _DefaultBuildParameters = GrainBoundaryWithVacuumBuilderParameters()
    _ConfigurationType: Type[GrainBoundaryWithVacuumConfiguration] = GrainBoundaryWithVacuumConfiguration
    _GeneratedItemType: Material = Material

    def _generate(self, configuration: GrainBoundaryWithVacuumConfiguration) -> List[Material]:
        # Create vacuum configuration
        vacuum_config = VacuumConfiguration(
            size=configuration.vacuum,
            crystal=configuration.grain_boundary_material,
            direction=AxisEnum.z,
        )

        # Stack grain boundary material with vacuum
        stack_config = StackConfiguration(
            stack_components=[configuration.grain_boundary_material, vacuum_config],
            direction=AxisEnum.z,
        )

        # Build the stack
        stack_builder = Stack2ComponentsBuilder()
        grain_boundary_with_vacuum = stack_builder.get_material(stack_config)

        # Apply supercell if specified
        if configuration.xy_supercell_matrix != [[1, 0], [0, 1]]:
            from mat3ra.made.tools.operations.core.unary import supercell

            grain_boundary_with_vacuum = supercell(grain_boundary_with_vacuum, configuration.xy_supercell_matrix)

        return [grain_boundary_with_vacuum]

    def _update_material_name(
        self, material: Material, configuration: GrainBoundaryWithVacuumConfiguration
    ) -> Material:
        formula = get_chemical_formula(configuration.grain_boundary_material)
        miller_indices_str = "".join([str(i) for i in configuration.miller_indices])
        new_name = f"{formula}({miller_indices_str}), Grain Boundary with Vacuum"
        material.name = new_name
        return material
