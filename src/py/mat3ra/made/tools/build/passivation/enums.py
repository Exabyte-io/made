# TODO: Get this from periodic table
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