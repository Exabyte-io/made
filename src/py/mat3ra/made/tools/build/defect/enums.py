from enum import Enum

from mat3ra.utils.factory import BaseFactory


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


class DefectBuilderFactory(BaseFactory):
    __class_registry__ = {
        "vacancy": "mat3ra.made.tools.build.defect.builders.VacancyPointDefectBuilder",
        "substitution": "mat3ra.made.tools.build.defect.builders.SubstitutionPointDefectBuilder",
        "interstitial": "mat3ra.made.tools.build.defect.builders.InterstitialPointDefectBuilder",
        "adatom:coordinate": "mat3ra.made.tools.build.defect.builders.AdatomSlabDefectBuilder",
        "adatom:crystal_site": "mat3ra.made.tools.build.defect.builders.CrystalSiteAdatomSlabDefectBuilder",
        "adatom:equidistant": "mat3ra.made.tools.build.defect.builders.EquidistantAdatomSlabDefectBuilder",
    }
