from enum import Enum


class PointDefectTypeEnum(str, Enum):
    VACANCY = "vacancy"
    SUBSTITUTION = "substitution"
    INTERSTITIAL = "interstitial"
    ADATOM = "adatom"


class ComplexDefectTypeEnum(str, Enum):
    PAIR = "pair"


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


class VacancyPlacementMethodEnum(Enum):
    CLOSEST_SITE = AtomPlacementMethodEnum.CLOSEST_SITE


class SubstitutionPlacementMethodEnum(Enum):
    CLOSEST_SITE = AtomPlacementMethodEnum.CLOSEST_SITE


class InterstitialPlacementMethodEnum(Enum):
    EXACT_COORDINATE = AtomPlacementMethodEnum.EXACT_COORDINATE
    VORONOI_SITE = AtomPlacementMethodEnum.VORONOI_SITE


class AdatomPlacementMethodEnum(Enum):
    EXACT_COORDINATE = AtomPlacementMethodEnum.EXACT_COORDINATE
    EQUIDISTANT = AtomPlacementMethodEnum.EQUIDISTANT
    NEW_CRYSTAL_SITE = AtomPlacementMethodEnum.NEW_CRYSTAL_SITE


class CoordinatesShapeEnum(str, Enum):
    SPHERE = "sphere"
    CYLINDER = "cylinder"
    BOX = "rectangle"
    TRIANGULAR_PRISM = "triangular_prism"
