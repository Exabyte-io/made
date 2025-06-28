from typing import Tuple, Optional, Union

from mat3ra.made.material import Material
from mat3ra.made.tools.build.slab.entities import Termination
from .builders import NanoribbonBuilder
from .configuration import (
    NanoribbonConfiguration, 
    get_miller_indices_from_edge_type,
    create_nanoribbon_configuration_from_edge_type
)


def create_nanoribbon(
    material: Material,
    miller_indices_uv: Union[Tuple[int, int], str],
    width: int,
    length: int,
    vacuum_width: int = 3,
    vacuum_length: int = 0,
    vacuum_z: float = 10.0,
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
        vacuum_width: The width of the vacuum region in number of unit cells.
        vacuum_length: The length of the vacuum region in number of unit cells.
        vacuum_z: The vacuum in the z direction in Angstroms.
        termination: The termination to use for the nanoribbon. If None, uses default termination.

    Returns:
        Material: The generated nanoribbon material.

    Example:
        # Using Miller indices
        nanoribbon = create_nanoribbon(monolayer, (1, 1), 10, 5, 3, 0, 10.0)
        
        # Using edge type shorthand
        zigzag_nanoribbon = create_nanoribbon(monolayer, "zigzag", 10, 5, 3, 0, 10.0)
        armchair_nanoribbon = create_nanoribbon(monolayer, "armchair", 10, 5, 3, 0, 10.0)
    """
    
    # Handle string input for edge types
    if isinstance(miller_indices_uv, str):
        miller_indices_uv = get_miller_indices_from_edge_type(miller_indices_uv)
    
    # Create configuration
    configuration = NanoribbonConfiguration.from_parameters(
        monolayer=material,
        miller_indices_uv=miller_indices_uv,
        width=width,
        length=length,
        termination=termination,
        vacuum_width=vacuum_width,
        vacuum_length=vacuum_length,
        vacuum_z=vacuum_z,
    )
    
    # Build the nanoribbon
    builder = NanoribbonBuilder()
    nanoribbon = builder.get_material(configuration)
    
    return nanoribbon


def create_nanoribbon_from_edge_type(
    material: Material,
    edge_type: str,
    width: int,
    length: int,
    vacuum_width: int = 3,
    vacuum_length: int = 0,
    vacuum_z: float = 10.0,
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
        vacuum_width: The width of the vacuum region in number of unit cells.
        vacuum_length: The length of the vacuum region in number of unit cells.
        vacuum_z: The vacuum in the z direction in Angstroms.
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
        vacuum_z=vacuum_z,
        termination=termination,
    )


def get_nanoribbon_terminations(
    material: Material, 
    miller_indices_uv: Union[Tuple[int, int], str]
) -> list:
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