from typing import List, Union, Optional, Tuple
from mat3ra.made.material import Material
from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum

from mat3ra.made.tools.build import BaseConfigurationPydantic
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


class NanoribbonConfiguration(StackConfiguration):
    """
    Configuration for building a nanoribbon from a monolayer material.

    Args:
        stack_components: List of configuration objects for nanoribbon components.
        direction: Direction along which to stack components.
        miller_indices_uv: The (u,v) Miller indices for the nanoribbon direction.
        width: The width of the nanoribbon in number of unit cells.
        length: The length of the nanoribbon in number of unit cells.
        vacuum_width: The width of the vacuum region in number of unit cells.
        vacuum_length: The length of the vacuum region in number of unit cells.
        vacuum_z: The vacuum in the z direction in Angstroms.
    """

    type: str = "NanoribbonConfiguration"
    stack_components: List[Union[CrystalLatticeLinesUniqueRepeatedConfiguration, VacuumConfiguration]]
    direction: AxisEnum = AxisEnum.x

    # Nanoribbon specific parameters
    miller_indices_uv: Tuple[int, int]
    width: int
    length: int
    vacuum_width: int = 3
    vacuum_length: int = 0
    vacuum_z: float = 10.0

    @property
    def lattice_lines(self):
        """Get the lattice lines configuration component."""
        return self.stack_components[0]

    @property
    def vacuum_configuration(self) -> VacuumConfiguration:
        """Get the vacuum configuration component."""
        return self.stack_components[1]

    @classmethod
    def from_parameters(
        cls,
        monolayer: Material,
        miller_indices_uv: Tuple[int, int],
        width: int,
        length: int,
        termination: Optional[Termination] = None,
        vacuum_width: int = 3,
        vacuum_length: int = 0,
        vacuum_z: float = 10.0,
    ) -> "NanoribbonConfiguration":
        """
        Create a nanoribbon configuration from parameters.

        Args:
            monolayer: The monolayer material to create the nanoribbon from.
            miller_indices_uv: The (u,v) Miller indices for the nanoribbon direction.
            width: The width of the nanoribbon in number of unit cells.
            length: The length of the nanoribbon in number of unit cells.
            termination: The termination to use for the nanoribbon.
            vacuum_width: The width of the vacuum region in number of unit cells.
            vacuum_length: The length of the vacuum region in number of unit cells.
            vacuum_z: The vacuum in the z direction in Angstroms.
        """
        from ...analyze.lattice_lines import CrystalLatticeLinesAnalyzer
        from ..nanoribbon.builders import CrystalLatticeLinesRepeatedBuilder

        # Get terminations if not provided
        if termination is None:
            analyzer = CrystalLatticeLinesAnalyzer(monolayer, miller_indices_uv)
            termination = analyzer.default_termination

        # Create lattice lines configuration
        lattice_lines_config = CrystalLatticeLinesUniqueRepeatedConfiguration(
            crystal=monolayer,
            miller_indices_uv=miller_indices_uv,
            termination_top=termination,
            number_of_repetitions_width=width,
            number_of_repetitions_length=length,
        )

        # Build the lattice lines material using the dedicated builder for vacuum configuration
        lattice_lines_builder = CrystalLatticeLinesRepeatedBuilder()
        lattice_lines_material = lattice_lines_builder.get_material(lattice_lines_config)

        # Create vacuum configuration
        vacuum_config = VacuumConfiguration(
            size=vacuum_z,
            crystal=lattice_lines_material,  # Use the built material for proper vacuum calculation
            direction=AxisEnum.z,
        )

        return cls(
            stack_components=[lattice_lines_config, vacuum_config],
            direction=AxisEnum.z,
            miller_indices_uv=miller_indices_uv,
            width=width,
            length=length,
            vacuum_width=vacuum_width,
            vacuum_length=vacuum_length,
            vacuum_z=vacuum_z,
        )


# Helper functions for shorthand notation
def get_miller_indices_from_edge_type(edge_type: str) -> Tuple[int, int]:
    """
    Convert edge type shorthand to (u,v) Miller indices.

    Args:
        edge_type: "zigzag" or "armchair"

    Returns:
        Tuple of (u,v) Miller indices.
    """
    if edge_type.lower() == "zigzag":
        return (1, 1)
    elif edge_type.lower() == "armchair":
        return (0, 1)
    else:
        raise ValueError(f"Unknown edge type: {edge_type}. Use 'zigzag' or 'armchair'.")


def create_nanoribbon_configuration_from_edge_type(
    monolayer: Material,
    edge_type: str,
    width: int,
    length: int,
    vacuum_width: int = 3,
    vacuum_length: int = 0,
    vacuum_z: float = 10.0,
) -> NanoribbonConfiguration:
    """
    Create a nanoribbon configuration using edge type shorthand.

    Args:
        monolayer: The monolayer material to create the nanoribbon from.
        edge_type: "zigzag" or "armchair"
        width: The width of the nanoribbon in number of unit cells.
        length: The length of the nanoribbon in number of unit cells.
        vacuum_width: The width of the vacuum region in number of unit cells.
        vacuum_length: The length of the vacuum region in number of unit cells.
        vacuum_z: The vacuum in the z direction in Angstroms.
    """
    miller_indices_uv = get_miller_indices_from_edge_type(edge_type)
    return NanoribbonConfiguration.from_parameters(
        monolayer=monolayer,
        miller_indices_uv=miller_indices_uv,
        width=width,
        length=length,
        vacuum_width=vacuum_width,
        vacuum_length=vacuum_length,
        vacuum_z=vacuum_z,
    )
