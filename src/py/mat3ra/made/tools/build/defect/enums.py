from enum import Enum


class PointDefectTypeEnum(str, Enum):
    VACANCY = "vacancy"
    SUBSTITUTION = "substitution"
    INTERSTITIAL = "interstitial"
    ADATOM = "adatom"
    PAIR = "pair"


class SlabDefectTypeEnum(str, Enum):
    ISLAND = "island"
    TERRACE = "terrace"


class AtomPlacementMethodEnum(str, Enum):
    COORDINATE = "coordinate"
    CLOSEST_SITE = "closest_site"
    EQUIDISTANT = "equidistant"
    CRYSTAL_SITE = "crystal_site"
