from enum import Enum


class EdgeTypes(str, Enum):
    """
    Enum for nanoribbon edge types.
    """

    zigzag = "zigzag"
    armchair = "armchair"
