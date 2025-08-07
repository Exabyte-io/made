from enum import Enum
from typing import Tuple


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
