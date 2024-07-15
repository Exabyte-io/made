from enum import Enum


class PointDefectTypeEnum(str, Enum):
    VACANCY = "vacancy"
    SUBSTITUTION = "substitution"
    INTERSTITIAL = "interstitial"
    ADATOM = "adatom"


class SlabDefectTypeEnum(str, Enum):
    ISLAND = "island"
