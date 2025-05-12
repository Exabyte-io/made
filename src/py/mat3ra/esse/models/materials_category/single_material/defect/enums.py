from enum import Enum


class PointDefectTypeEnum(str, Enum):
    """Enum for point defect types."""
    VACANCY = "vacancy"
    SUBSTITUTION = "substitution"
    INTERSTITIAL = "interstitial"
    ADATOM = "adatom"


class AtomPlacementMethodEnum(str, Enum):
    """Enum for atom placement methods."""
    COORDINATE = "coordinate"
    VORONOI_SITE = "voronoi_site" 