from typing import Type, Optional

from mat3ra.made.material import Material
from .configuration import GrainBoundaryLinearConfiguration
from .configuration import GrainBoundaryPlanarConfiguration
from ..interface.builders import InterfaceBuilder, InterfaceBuilderParameters
from ...modify import wrap_to_unit_cell
from ...operations.core.unary import supercell


class GrainBoundaryBuilderParameters(InterfaceBuilderParameters):
    pass


class GrainBoundaryLinearBuilderParameters(GrainBoundaryBuilderParameters):
    max_supercell_matrix_int: Optional[int] = None
    limit_max_int: Optional[int] = 30
    angle_tolerance: float = 0.1
    return_first_match: bool = False
    edge_inclusion_tolerance: float = 1.0
    distance_tolerance: float = 1.0


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


class GrainBoundaryLinearBuilder(InterfaceBuilder):
    """
    Creates a linear grain boundary by stacking two materials along x or y direction.
    """

    _ConfigurationType = GrainBoundaryLinearConfiguration
    _BuildParametersType = GrainBoundaryLinearBuilderParameters
    _DefaultBuildParameters = GrainBoundaryLinearBuilderParameters()

    def get_name_suffix(self, configuration: GrainBoundaryLinearConfiguration) -> str:
        angle_str = f"{configuration.actual_angle:.2f} degrees" if configuration.actual_angle is not None else ""
        return f"Linear Grain Boundary, {angle_str}"
