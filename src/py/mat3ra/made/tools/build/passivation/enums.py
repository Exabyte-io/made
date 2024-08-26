# TODO: Get this from peroidic table
from enum import Enum

BOND_LENGTHS_MAP = {
    ("C", "H"): 1.09,
    ("Ni", "H"): 1.09,
    ("Si", "H"): 1.48,
}


class SurfaceTypes(str, Enum):
    TOP = "top"
    BOTTOM = "bottom"
    BOTH = "both"


class EdgeTypes(str, Enum):
    ALONG_A = "along_a"
    ALONG_B = "along_b"
    BOTH = "both"
