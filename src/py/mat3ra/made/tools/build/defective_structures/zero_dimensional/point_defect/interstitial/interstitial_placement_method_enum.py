from enum import Enum

from ..atom_placement_method_enum import AtomPlacementMethodEnum


class InterstitialPlacementMethodEnum(Enum):
    EXACT_COORDINATE = AtomPlacementMethodEnum.EXACT_COORDINATE
    VORONOI_SITE = AtomPlacementMethodEnum.VORONOI_SITE
