from enum import Enum


class PointDefectTypeEnum(str, Enum):
    VACANCY = "vacancy"
    SUBSTITUTION = "substitution"
    INTERSTITIAL = "interstitial"
    ADATOM = "adatom"


class SlabDefectTypeEnum(str, Enum):
    ISLAND = "island"
    TERRACE = "terrace"


class ComplexDefectTypeEnum(str, Enum):
    PAIR = "pair"


class AtomPlacementMethodEnum(str, Enum):
    # Places the atom at the exact given coordinate.
    COORDINATE = "coordinate"
    # Among existing atoms, selects the closest one to the given coordinate.
    CLOSEST_SITE = "closest_site"
    # Places the atom at the equal distance from the closest atoms to the given coordinate.
    EQUIDISTANT = "equidistant"
    # Places the atom at the existing or extrapolated crystal site closest to the given coordinate.
    CRYSTAL_SITE = "crystal_site"


class CoordinatesShapeEnum(str, Enum):
    SPHERE = "sphere"
    CYLINDER = "cylinder"
    BOX = "rectangle"
    TRIANGULAR_PRISM = "triangular_prism"
