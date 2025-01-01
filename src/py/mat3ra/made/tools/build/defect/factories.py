from mat3ra.utils.factory import BaseFactory


class DefectBuilderFactory(BaseFactory):
    __class_registry__ = {
        "vacancy": "mat3ra.made.tools.build.defect.builders.VacancyPointDefectBuilder",
        "substitution": "mat3ra.made.tools.build.defect.builders.SubstitutionPointDefectBuilder",
        "interstitial": "mat3ra.made.tools.build.defect.builders.InterstitialPointDefectBuilder",
        "interstitial:voronoi_site": "mat3ra.made.tools.build.defect.builders.VoronoiInterstitialPointDefectBuilder",
        "adatom:coordinate": "mat3ra.made.tools.build.defect.builders.AdatomSlabDefectBuilder",
        "adatom:crystal_site": "mat3ra.made.tools.build.defect.builders.CrystalSiteAdatomSlabDefectBuilder",
        "adatom:equidistant": "mat3ra.made.tools.build.defect.builders.EquidistantAdatomSlabDefectBuilder",
    }
