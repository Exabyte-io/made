from typing import Type

from mat3ra.made.material import Material
from mat3ra.made.tools.build.interface.builders import InterfaceBuilder, InterfaceBuilderParameters
from .configuration import GrainBoundaryConfiguration
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
        # Since GrainBoundaryConfiguration now inherits from InterfaceConfiguration,
        # we can use the parent's _generate method directly
        return super()._generate(configuration)

    def _update_material_name(self, material: Material, configuration: GrainBoundaryConfiguration) -> Material:
        phase_1_formula = get_chemical_formula(configuration.phase_1_configuration.atomic_layers.crystal)
        phase_2_formula = get_chemical_formula(configuration.phase_2_configuration.atomic_layers.crystal)
        phase_1_miller = "".join([str(i) for i in configuration.phase_1_configuration.atomic_layers.miller_indices])
        phase_2_miller = "".join([str(i) for i in configuration.phase_2_configuration.atomic_layers.miller_indices])

        new_name = f"{phase_1_formula}({phase_1_miller})-{phase_2_formula}({phase_2_miller}), Grain Boundary"
        material.name = new_name
        return material
