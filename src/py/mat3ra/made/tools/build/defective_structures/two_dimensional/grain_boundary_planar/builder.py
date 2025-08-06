from typing import Type

from mat3ra.made.tools.build.grain_boundary.builders import GrainBoundaryBuilderParameters

from mat3ra.made.material import Material
from mat3ra.made.tools.build.compound_pristine_structures.two_dimensional.interface import InterfaceBuilder
from mat3ra.made.tools.build.defective_structures.two_dimensional.grain_boundary_planar.configuration import (
    GrainBoundaryPlanarConfiguration,
)
from mat3ra.made.tools.modify import wrap_to_unit_cell
from mat3ra.made.tools.operations.core.unary import supercell


class GrainBoundaryPlanarBuilder(InterfaceBuilder):
    """
    Builder for creating grain boundaries.

    Uses InterfaceBuilder to create a Z-stacked interface, then transforms the supercell
    to flip the Z direction to X direction for grain boundary orientation.
    """

    _BuildParametersType: Type[GrainBoundaryBuilderParameters] = GrainBoundaryBuilderParameters
    _DefaultBuildParameters = GrainBoundaryBuilderParameters()

    def _generate(self, configuration: GrainBoundaryPlanarConfiguration) -> Material:
        interface = super()._generate(configuration)
        rotated_interface = supercell(interface, [[0, 0, 1], [0, 1, 0], [1, 0, 0]])
        wrapped_interface = wrap_to_unit_cell(rotated_interface)
        return wrapped_interface

    def get_name_suffix(self, configuration: GrainBoundaryPlanarConfiguration) -> str:
        return "Grain Boundary"
