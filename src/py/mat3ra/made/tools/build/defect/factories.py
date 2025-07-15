from typing import List, Dict, Union, Type

from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.crystal_site import CrystalSiteAnalyzer, VoronoiCrystalSiteAnalyzer
from mat3ra.made.tools.build.defect.enums import (
    PointDefectTypeEnum,
    VacancyPlacementMethodEnum,
    SubstitutionPlacementMethodEnum,
    InterstitialPlacementMethodEnum,
    AtomPlacementMethodEnum,
)
from mat3ra.made.tools.build.defect.point.configuration import (
    VacancyDefectConfiguration,
    SubstitutionalDefectConfiguration,
    InterstitialDefectConfiguration,
    PointDefectConfiguration,
)


def resolve_coordinate(material: Material, coordinate: List[float], placement_method) -> List[float]:
    if placement_method in [VacancyPlacementMethodEnum.CLOSEST_SITE, SubstitutionPlacementMethodEnum.CLOSEST_SITE]:
        analyzer = CrystalSiteAnalyzer(material=material, coordinate=coordinate)
        return analyzer.closest_site_coordinate
    elif placement_method == InterstitialPlacementMethodEnum.VORONOI_SITE:
        analyzer = VoronoiCrystalSiteAnalyzer(material=material, coordinate=coordinate)
        return analyzer.voronoi_site_coordinate
    else:
        analyzer = CrystalSiteAnalyzer(material=material, coordinate=coordinate)
        return analyzer.exact_coordinate


DEFECT_CONFIGURATION_FACTORIES: Dict[PointDefectTypeEnum, Type[PointDefectConfiguration]] = {
    PointDefectTypeEnum.VACANCY: VacancyDefectConfiguration,
    PointDefectTypeEnum.SUBSTITUTION: SubstitutionalDefectConfiguration,
    PointDefectTypeEnum.INTERSTITIAL: InterstitialDefectConfiguration,
}


def create_defect_configuration(
    material: Material,
    defect_type: PointDefectTypeEnum,
    coordinate: List[float],
    element: str = None,
    placement_method: Union[
        VacancyPlacementMethodEnum, SubstitutionPlacementMethodEnum, InterstitialPlacementMethodEnum
    ] = AtomPlacementMethodEnum.EXACT_COORDINATE,
) -> PointDefectConfiguration:
    if defect_type not in DEFECT_CONFIGURATION_FACTORIES:
        raise ValueError(f"Unknown defect type: {defect_type}")

    resolved_coordinate = resolve_coordinate(material, coordinate, placement_method)

    config_class = DEFECT_CONFIGURATION_FACTORIES[defect_type]
    return config_class.from_parameters(crystal=material, coordinate=resolved_coordinate, element=element)
