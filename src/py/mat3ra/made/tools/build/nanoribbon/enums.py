from enum import Enum


class EdgeTypes(str, Enum):
    """
    Enum for edge types.
    """

    zigzag = "zigzag"
    armchair = "armchair"
