from typing import Tuple, Optional, Union

from mat3ra.made.material import Material
from mat3ra.made.tools.build.slab.entities import Termination
from .builders import NanoribbonBuilder, NanoTapeBuilder, NanoTapeBuilderParameters, NanoribbonBuilderParameters
from .configuration import (
    NanoribbonConfiguration,
    get_miller_indices_from_edge_type,
)
from mat3ra.made.tools.analyze.nanoribbon_analyzer import NanoribbonAnalyzer
from mat3ra.made.tools.analyze.nanotape_analyzer import NanoTapeAnalyzer


def create_nanotape(
    material: Material,
    miller_indices_uv: Union[Tuple[int, int], str],
    width: int,
    length: int,
    vacuum_width: float = 10.0,
    use_rectangular_cell: bool = True,
) -> Material:
    """
    Create a nanotape from a monolayer material.

    Args:
        material: The monolayer material to create the nanotape from.
        miller_indices_uv: The (u,v) Miller indices or edge type ("zigzag"/"armchair").
        width: The width of the nanotape in number of unit cells.
        length: The length of the nanotape in number of unit cells.
        vacuum_width: The width of the vacuum region in Angstroms (cartesian).
        use_rectangular_cell: Whether the nanotape is rectangular.

    Returns:
        The nanotape material.
    """
    # Convert edge type to Miller indices if needed
    if isinstance(miller_indices_uv, str):
        miller_indices_uv = get_miller_indices_from_edge_type(miller_indices_uv)

    # Create configuration using analyzer
    analyzer = NanoTapeAnalyzer(material, miller_indices_uv)
    configuration = analyzer.get_configuration(
        width=width,
        length=length,
        vacuum_width=vacuum_width,
    )

    # Build material using builder
    build_parameters = NanoTapeBuilderParameters(use_rectangular_lattice=use_rectangular_cell)
    builder = NanoTapeBuilder(
        build_parameters=build_parameters,
    )
    return builder.get_material(configuration)


def create_nanoribbon(
    material: Material,
    miller_indices_uv: Union[Tuple[int, int], str],
    width: int,
    length: int,
    vacuum_width: float = 10.0,
    vacuum_length: float = 0.0,
    use_rectangular_cell: bool = True,
    termination: Optional[Termination] = None,
) -> Material:
    """
    Creates a nanoribbon material from a monolayer material.

    Args:
        material: The monolayer material to create the nanoribbon from (assumes vacuum is present).
        miller_indices_uv: The (u,v) Miller indices for the nanoribbon direction, or edge type string
                          ("zigzag" for (1,1), "armchair" for (0,1)).
        width: The width of the nanoribbon in number of unit cells.
        length: The length of the nanoribbon in number of unit cells.
        vacuum_width: The width of the vacuum region in Angstroms (cartesian).
        vacuum_length: The length of the vacuum region in Angstroms (cartesian).
        use_rectangular_cell: Whether the nanoribbon is rectangular.
        termination: The termination to use for the nanoribbon. If None, uses default termination.

    Returns:
        Material: The generated nanoribbon material.

    Example:
        # Using Miller indices
        nanoribbon = create_nanoribbon(monolayer, (1, 1), 10, 5, 10.0, 0.0)

        # Using edge type shorthand
        zigzag_nanoribbon = create_nanoribbon(monolayer, "zigzag", 10, 5, 10.0, 0.0)
        armchair_nanoribbon = create_nanoribbon(monolayer, "armchair", 10, 5, 10.0, 0.0)
    """

    # Handle string input for edge types
    if isinstance(miller_indices_uv, str):
        miller_indices_uv = get_miller_indices_from_edge_type(miller_indices_uv)

    # Create configuration using analyzer
    analyzer = NanoribbonAnalyzer(material, miller_indices_uv)
    configuration = analyzer.get_configuration(
        width=width,
        length=length,
        vacuum_width=vacuum_width,
        vacuum_length=vacuum_length,
    )

    # Build material using builder
    build_parameters = NanoribbonBuilderParameters(use_rectangular_lattice=use_rectangular_cell)
    builder = NanoribbonBuilder(
        build_parameters=build_parameters,
    )
    nanoribbon = builder.get_material(configuration)

    return nanoribbon


def create_nanoribbon_from_edge_type(
    material: Material,
    edge_type: str,
    width: int,
    length: int,
    vacuum_width: float = 10.0,
    vacuum_length: float = 0.0,
    termination: Optional[Termination] = None,
) -> Material:
    """
    Creates a nanoribbon material using edge type shorthand.

    This is a convenience function that wraps create_nanoribbon for edge type strings.

    Args:
        material: The monolayer material to create the nanoribbon from.
        edge_type: "zigzag" or "armchair"
        width: The width of the nanoribbon in number of unit cells.
        length: The length of the nanoribbon in number of unit cells.
        vacuum_width: The width of the vacuum region in Angstroms (cartesian).
        vacuum_length: The length of the vacuum region in Angstroms (cartesian).
        termination: The termination to use for the nanoribbon.

    Returns:
        Material: The generated nanoribbon material.
    """
    return create_nanoribbon(
        material=material,
        miller_indices_uv=edge_type,
        width=width,
        length=length,
        vacuum_width=vacuum_width,
        vacuum_length=vacuum_length,
        termination=termination,
    )


def get_nanoribbon_terminations(material: Material, miller_indices_uv: Union[Tuple[int, int], str]) -> list:
    """
    Get available terminations for a nanoribbon from a monolayer material.

    Args:
        material: The monolayer material.
        miller_indices_uv: The (u,v) Miller indices or edge type string.

    Returns:
        List of available terminations.
    """
    from mat3ra.made.tools.analyze.lattice_lines import CrystalLatticeLinesAnalyzer

    # Handle string input for edge types
    if isinstance(miller_indices_uv, str):
        miller_indices_uv = get_miller_indices_from_edge_type(miller_indices_uv)

    analyzer = CrystalLatticeLinesAnalyzer(material, miller_indices_uv)
    return analyzer.terminations
