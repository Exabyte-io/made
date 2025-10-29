from enum import Enum
from typing import Tuple


class EdgeTypesEnum(str, Enum):
    """
    Enum for nanoribbon/nanotape edge types.
    """

    zigzag = "zigzag"
    armchair = "armchair"


MILLER_INDICES_TO_EDGE_TYPE_MAP = {
    (0, 1): EdgeTypesEnum.zigzag,
    (1, 1): EdgeTypesEnum.armchair,
}


def get_miller_indices_from_edge_type(edge_type: EdgeTypesEnum) -> Tuple[int, int]:
    """
    Convert edge type shorthand to (u,v) Miller indices.

    Args:
        edge_type: "zigzag" or "armchair"

    Returns:
        Tuple of (u,v) Miller indices.
    """
    for miller_indices, mapped_edge_type in MILLER_INDICES_TO_EDGE_TYPE_MAP.items():
        if mapped_edge_type == edge_type:
            return miller_indices
    raise ValueError(f"Unknown edge type: {edge_type}. Use 'zigzag' or 'armchair'.")


def get_edge_type_from_miller_indices(miller_indices_2d: Tuple[int, int]) -> str:
    """
    Convert (u,v) Miller indices to edge type string.

    Args:
        miller_indices_2d: Tuple of (u,v) Miller indices

    Returns:
        Edge type string: "Zigzag", "Armchair", or "(uv)" for other indices
    """
    edge_type = MILLER_INDICES_TO_EDGE_TYPE_MAP.get(miller_indices_2d)
    if edge_type:
        return edge_type.value.capitalize()
    else:
        miller_str = f"{miller_indices_2d[0]}{miller_indices_2d[1]}"
        return f"({miller_str})"
