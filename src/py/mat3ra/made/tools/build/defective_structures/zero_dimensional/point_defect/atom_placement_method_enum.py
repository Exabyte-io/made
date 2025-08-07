from enum import Enum


class AtomPlacementMethodEnum(str, Enum):
    # Places the atom at the exact given coordinate.
    EXACT_COORDINATE = "exact_coordinate"
    # Among existing atoms, selects the closest one to the given coordinate.
    CLOSEST_SITE = "closest_site"
    # Places the atom at the equal distance from the closest atoms to the given coordinate.
    EQUIDISTANT = "equidistant"
    # Places the atom at the existing or extrapolated crystal site closest to the given coordinate.
    NEW_CRYSTAL_SITE = "new_crystal_site"
    # Places the atom at Voronoi site closest to the given coordinate.
    VORONOI_SITE = "voronoi_site"
