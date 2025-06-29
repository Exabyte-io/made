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

    @classmethod
    def from_parameters(
        cls,
        material: Material,
        miller_indices_uv: Optional[Tuple[int, int]] = None,
        edge_type: Optional[EdgeTypes] = EdgeTypes.zigzag,
        width: int = 2,
        length: int = 2,
        vacuum_width: float = 10.0,
        termination: Optional[Termination] = None,
    ) -> "NanoTapeConfiguration":
        """
        Create a NanoTapeConfiguration from parameters.

        Args:
            material: The monolayer material to create the nanotape from.
            miller_indices_uv: The (u,v) Miller indices for the nanotape direction.
            edge_type: Edge type string ("zigzag"/"armchair"). Optional if miller_indices_uv is provided.
            width: The width of the nanotape in number of unit cells.
            length: The length of the nanotape in number of unit cells.
            vacuum_width: The width of the vacuum region in Angstroms (cartesian).
            termination: The termination to use for the nanotape. If None, uses default termination.

        Returns:
            NanoTapeConfiguration: The created configuration.
        """
        from mat3ra.made.tools.analyze.nanotape_analyzer import NanoTapeAnalyzer

        if miller_indices_uv is None and edge_type is None:
            raise ValueError("Either miller_indices_uv or edge_type must be provided")

        if miller_indices_uv is None and edge_type is not None:
            miller_indices_uv = get_miller_indices_from_edge_type(edge_type)

        analyzer = NanoTapeAnalyzer(material=material, miller_indices_uv=miller_indices_uv)
        return analyzer.get_configuration(
            width=width,
            length=length,
            vacuum_width=vacuum_width,
            termination=termination,
        )


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

    @classmethod
    def from_parameters(
        cls,
        material: Material,
        miller_indices_uv: Optional[Tuple[int, int]] = None,
        edge_type: Optional[EdgeTypes] = EdgeTypes.zigzag,
        width: int = 2,
        length: int = 2,
        vacuum_width: float = 10.0,
        vacuum_length: float = 10.0,
        termination: Optional[Termination] = None,
    ) -> "NanoribbonConfiguration":
        """
        Create a NanoribbonConfiguration from parameters.

        Args:
            material: The monolayer material to create the nanoribbon from.
            miller_indices_uv: The (u,v) Miller indices for the nanoribbon direction.
            edge_type: Edge type string ("zigzag"/"armchair"). Optional if miller_indices_uv is provided.
            width: The width of the nanoribbon in number of unit cells.
            length: The length of the nanoribbon in number of unit cells.
            vacuum_width: The width of the vacuum region in Angstroms (cartesian).
            vacuum_length: The length of the vacuum region in Angstroms (cartesian).
            termination: The termination to use for the nanoribbon. If None, uses default termination.

        Returns:
            NanoribbonConfiguration: The created configuration.
        """
        from mat3ra.made.tools.analyze.nanoribbon_analyzer import NanoribbonAnalyzer

        if miller_indices_uv is None and edge_type is None:
            raise ValueError("Either miller_indices_uv or edge_type must be provided")

        if miller_indices_uv is None and edge_type is not None:
            miller_indices_uv = get_miller_indices_from_edge_type(edge_type)

        analyzer = NanoribbonAnalyzer(material=material, miller_indices_uv=miller_indices_uv)
        return analyzer.get_configuration(
            width=width,
            length=length,
            vacuum_width=vacuum_width,
            vacuum_length=vacuum_length,
            termination=termination,
        )


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
