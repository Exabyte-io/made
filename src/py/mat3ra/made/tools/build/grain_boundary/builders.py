from typing import Type

from mat3ra.made.material import Material
from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from ..interface.builders import InterfaceBuilder, InterfaceBuilderParameters
from ...convert.utils import InterfacePartsEnum
from .configuration import GrainBoundaryConfiguration
from ...analyze.other import get_chemical_formula
from ...operations.core.unary import supercell
from ...modify import wrap_to_unit_cell


class GrainBoundaryBuilderParameters(InterfaceBuilderParameters):
    pass


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
