from typing import Tuple, Optional, Union
from mat3ra.made.material import Material
from mat3ra.made.tools.build import BaseConfigurationPydantic, MaterialWithBuildMetadata
from ..slab.entities import Termination
from enum import Enum


class EdgeTypes(str, Enum):
    """
    Enum for nanoribbon/nanotape edge types.
    """

    zigzag = "zigzag"
    armchair = "armchair"


def get_miller_indices_from_edge_type(edge_type: EdgeTypes) -> Tuple[int, int]:
    """
    Convert edge type shorthand to (u,v) Miller indices.

    Args:
        edge_type: "zigzag" or "armchair"

    Returns:
        Tuple of (u,v) Miller indices.
    """
    if edge_type == EdgeTypes.zigzag:
        return (1, 1)
    elif edge_type == EdgeTypes.armchair:
        return (0, 1)
    else:
        raise ValueError(f"Unknown edge type: {edge_type}. Use 'zigzag' or 'armchair'.")


class CrystalLatticeLinesConfiguration(BaseConfigurationPydantic):
    """
    Configuration for creating crystal lattice lines from a material.

    Args:
        crystal: The monolayer material to create the lattice lines from.
        miller_indices_2d: The (u,v) Miller indices for the line direction.
    """

    crystal: Union[Material, MaterialWithBuildMetadata]
    miller_indices_2d: Tuple[int, int]


class CrystalLatticeLinesUniqueRepeatedConfiguration(CrystalLatticeLinesConfiguration):
    """
    Configuration for creating repeated crystal lattice lines with termination.

    Args:
        crystal: The monolayer material to create the lattice lines from.
        miller_indices_2d: The (u,v) Miller indices for the line direction.
        termination_top: The termination to use for the lattice lines.
    """

    # TODO: right and left for x terminations
    termination_top: Termination
    termination_bottom: Optional[Termination] = None
    number_of_repetitions_width: int = 1
    number_of_repetitions_length: int = 1
