from typing import List, Optional

import numpy as np
from mat3ra.code.entity import InMemoryEntityPydantic
from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.interface.grain_boundary import GrainBoundaryAnalyzer
from mat3ra.made.tools.analyze.lattice_planes import CrystalLatticePlanesMaterialAnalyzer
from mat3ra.made.tools.build.slab.configurations import SlabConfiguration
from mat3ra.made.tools.build.slab.builders import SlabBuilder, SlabBuilderParameters
from mat3ra.made.tools.build.stack.builders import Stack2ComponentsBuilder
from mat3ra.made.tools.build.stack.configuration import StackConfiguration
from mat3ra.made.tools.build.slab.configurations import SlabStrainedSupercellWithGapConfiguration
from mat3ra.made.tools.build.vacuum.configuration import VacuumConfiguration
from mat3ra.made.tools.build import BaseBuilder, BaseBuilderParameters
from mat3ra.made.tools.build.utils import stack_two_materials_xy
from mat3ra.made.tools.modify import translate_by_vector
from ...analyze.other import get_chemical_formula

from .configuration import GrainBoundaryConfiguration, GrainBoundaryWithVacuumConfiguration


class GrainBoundaryBuilderParameters(BaseBuilderParameters):
    """Parameters for grain boundary builder."""
    pass


class GrainBoundaryBuilder(BaseBuilder):
    """
    Builder for creating grain boundaries.

    Creates a grain boundary by stacking two phases in x/y direction
    with a relative shift perpendicular to the interface.
    """

    _BuildParametersType: type(GrainBoundaryBuilderParameters) = GrainBoundaryBuilderParameters
    _DefaultBuildParameters = GrainBoundaryBuilderParameters()
    _ConfigurationType: type(GrainBoundaryConfiguration) = GrainBoundaryConfiguration
    _GeneratedItemType: Material = Material

    def _generate(self, configuration: GrainBoundaryConfiguration) -> List[Material]:
        # Get the strained configurations from the analyzer
        phase_1_config = configuration.phase_1_configuration
        phase_2_config = configuration.phase_2_configuration

        # Create materials from the strained configurations
        from mat3ra.made.tools.build.slab.builders import SlabStrainedSupercellBuilder

        phase_1_builder = SlabStrainedSupercellBuilder()
        phase_1_material = phase_1_builder.get_material(phase_1_config.substrate_configuration)

        phase_2_builder = SlabStrainedSupercellBuilder()
        phase_2_material = phase_2_builder.get_material(phase_2_config.film_configuration)

        # Apply translation to phase 2 if specified
        if any(v != 0.0 for v in configuration.translation_vector):
            phase_2_material = translate_by_vector(phase_2_material, configuration.translation_vector)

        # Stack the materials in x/y direction with gap
        grain_boundary = stack_two_materials_xy(
            phase_1_material,
            phase_2_material,
            gap=configuration.gap,
            edge_inclusion_tolerance=1.0,
            distance_tolerance=1.0,
        )

        return [grain_boundary]

    def _update_material_name(self, material: Material, configuration: GrainBoundaryConfiguration) -> Material:
        phase_1_formula = get_chemical_formula(configuration.phase_1_configuration.substrate_configuration.atomic_layers.crystal)
        phase_2_formula = get_chemical_formula(configuration.phase_2_configuration.film_configuration.atomic_layers.crystal)
        phase_1_miller = "".join([str(i) for i in configuration.phase_1_configuration.substrate_configuration.atomic_layers.miller_indices])
        phase_2_miller = "".join([str(i) for i in configuration.phase_2_configuration.film_configuration.atomic_layers.miller_indices])

        new_name = f"{phase_1_formula}({phase_1_miller})-{phase_2_formula}({phase_2_miller}), Grain Boundary"
        material.name = new_name
        return material


class GrainBoundaryWithVacuumBuilderParameters(BaseBuilderParameters):
    """Parameters for grain boundary with vacuum builder."""
    use_orthogonal_c: bool = True


class GrainBoundaryWithVacuumBuilder(BaseBuilder):
    """
    Builder for creating grain boundaries with vacuum (slab).

    Adds vacuum directly to the grain boundary material instead of trying to create
    a slab from it, since grain boundary materials don't have conventional terminations.
    """

    _BuildParametersType: type(GrainBoundaryWithVacuumBuilderParameters) = GrainBoundaryWithVacuumBuilderParameters
    _DefaultBuildParameters = GrainBoundaryWithVacuumBuilderParameters()
    _ConfigurationType: type(GrainBoundaryWithVacuumConfiguration) = GrainBoundaryWithVacuumConfiguration
    _GeneratedItemType: Material = Material

    def _generate(self, configuration: GrainBoundaryWithVacuumConfiguration) -> List[Material]:
        # Add vacuum directly to the grain boundary material
        from mat3ra.made.tools.build.vacuum.builders import VacuumBuilder
        from mat3ra.made.tools.build.stack.configuration import StackConfiguration
        from mat3ra.made.tools.build.stack.builders import Stack2ComponentsBuilder
        from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum

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

    def _update_material_name(self, material: Material, configuration: GrainBoundaryWithVacuumConfiguration) -> Material:
        formula = get_chemical_formula(configuration.grain_boundary_material)
        miller_indices_str = "".join([str(i) for i in configuration.miller_indices])
        new_name = f"{formula}({miller_indices_str}), Grain Boundary with Vacuum"
        material.name = new_name
        return material
