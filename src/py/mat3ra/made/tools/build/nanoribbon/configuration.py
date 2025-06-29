from typing import List, Union, Optional, Tuple
from mat3ra.made.material import Material
from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum

from mat3ra.made.tools.build import BaseConfigurationPydantic
from .enums import EdgeTypes
from ..stack.configuration import StackConfiguration
from ..vacuum.configuration import VacuumConfiguration
from ..slab.entities import Termination


class CrystalLatticeLinesConfiguration(BaseConfigurationPydantic):
    """
    Configuration for creating crystal lattice lines from a material.

    Args:
        crystal: The monolayer material to create the lattice lines from.
        miller_indices_uv: The (u,v) Miller indices for the line direction.
    """

    crystal: Material
    miller_indices_uv: Tuple[int, int]


class CrystalLatticeLinesUniqueRepeatedConfiguration(CrystalLatticeLinesConfiguration):
    """
    Configuration for creating repeated crystal lattice lines with termination.

    Args:
        crystal: The monolayer material to create the lattice lines from.
        miller_indices_uv: The (u,v) Miller indices for the line direction.
        termination_top: The termination to use for the lattice lines.
    """

    # TODO: right and left for x terminations
    termination_top: Termination
    termination_bottom: Optional[Termination] = None
    number_of_repetitions_width: int = 1
    number_of_repetitions_length: int = 1


class NanoTapeConfiguration(StackConfiguration):
    """
    Configuration for building a nanotape from crystal lattice lines.
    NanoTape = [CLLUR, vacuum] stacked on Y direction.

    Args:
        stack_components: List of configuration objects for nanotape components.
        direction: Direction along which to stack components.
    """

    type: str = "NanoTapeConfiguration"
    stack_components: List[Union[CrystalLatticeLinesUniqueRepeatedConfiguration, VacuumConfiguration]]
    direction: AxisEnum = AxisEnum.y
    use_rectangular_lattice: bool = True

    @property
    def lattice_lines(self):
        """Get the lattice lines configuration component."""
        return self.stack_components[0]

    @property
    def vacuum_configuration(self) -> VacuumConfiguration:
        """Get the vacuum configuration component."""
        return self.stack_components[1]


class NanoribbonConfiguration(StackConfiguration):
    """
    Configuration for building a nanoribbon from a nanotape.
    Nanoribbon = [NanoTape, vacuum] stacked on X direction.

    Args:
        stack_components: List of configuration objects for nanoribbon components.
        direction: Direction along which to stack components.
    """

    type: str = "NanoribbonConfiguration"
    stack_components: List[Union[NanoTapeConfiguration, VacuumConfiguration]]
    direction: AxisEnum = AxisEnum.x
    use_rectangular_lattice: bool = True

    @property
    def nanotape(self):
        """Get the nanotape configuration component."""
        return self.stack_components[0]

    @property
    def vacuum_configuration(self) -> VacuumConfiguration:
        """Get the vacuum configuration component."""
        return self.stack_components[1]


# Helper functions for shorthand notation
def get_miller_indices_from_edge_type(edge_type: EdgeTypes) -> Tuple[int, int]:
    """
    Convert edge type shorthand to (u,v) Miller indices.

    Args:
        edge_type: "zigzag" or "armchair"

    Returns:
        Tuple of (u,v) Miller indices.
    """
    if edge_type.lower() == EdgeTypes.zigzag:
        return (1, 1)
    elif edge_type.lower() == EdgeTypes.armchair:
        return (0, 1)
    else:
        raise ValueError(f"Unknown edge type: {edge_type}. Use 'zigzag' or 'armchair'.")
