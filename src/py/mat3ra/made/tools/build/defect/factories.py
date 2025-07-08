from mat3ra.utils.factory import BaseFactory


class DefectBuilderFactory(BaseFactory):
    __class_registry__ = {
        "adatom:coordinate": "mat3ra.made.tools.build.defect.builders.AdatomSlabDefectBuilder",
        "adatom:crystal_site": "mat3ra.made.tools.build.defect.builders.CrystalSiteAdatomSlabDefectBuilder",
        "adatom:equidistant": "mat3ra.made.tools.build.defect.builders.EquidistantAdatomSlabDefectBuilder",
    }
