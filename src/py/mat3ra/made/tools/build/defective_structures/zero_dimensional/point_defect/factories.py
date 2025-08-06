from typing import List, Dict, Union, Type, Optional

from mat3ra.made.material import Material
from .....analyze.crystal_site.crystal_site_analyzer import CrystalSiteAnalyzer
from .....analyze.crystal_site.voronoi_crystal_site_analyzer import VoronoiCrystalSiteAnalyzer
from .....build_components import MaterialWithBuildMetadata
from .atom_placement_method_enum import AtomPlacementMethodEnum
from .base.configuration import PointDefectConfiguration
from .interstitial.configuration import InterstitialDefectConfiguration
from .interstitial.interstitial_placement_method_enum import InterstitialPlacementMethodEnum
from .point_defect_type_enum import PointDefectTypeEnum
from .substitutional.configuration import SubstitutionalDefectConfiguration
from .substitutional.substitution_placement_method_enum import SubstitutionPlacementMethodEnum
from .vacancy.configuration import VacancyDefectConfiguration
from .vacancy.vacancy_placement_method_enum import VacancyPlacementMethodEnum


def resolve_coordinate(
    material: Union[Material, MaterialWithBuildMetadata],
    coordinate: List[float],
    placement_method,
    use_cartesian_coordinates: bool = False,
) -> List[float]:
    if use_cartesian_coordinates:
        coordinate = material.basis.cell.convert_point_to_crystal(coordinate)

    if placement_method in [VacancyPlacementMethodEnum.CLOSEST_SITE, SubstitutionPlacementMethodEnum.CLOSEST_SITE]:
        analyzer = CrystalSiteAnalyzer(material=material, coordinate=coordinate)
        return analyzer.closest_site_coordinate
    elif placement_method == InterstitialPlacementMethodEnum.VORONOI_SITE:
        analyzer = VoronoiCrystalSiteAnalyzer(material=material, coordinate=coordinate)
        return analyzer.voronoi_site_coordinate
    else:
        analyzer = CrystalSiteAnalyzer(material=material, coordinate=coordinate)
        return analyzer.exact_coordinate


class PointDefectConfigurationFactory:
    _type_to_configuration_map: Dict[PointDefectTypeEnum, Type[PointDefectConfiguration]] = {
        PointDefectTypeEnum.VACANCY: VacancyDefectConfiguration,
        PointDefectTypeEnum.SUBSTITUTION: SubstitutionalDefectConfiguration,
        PointDefectTypeEnum.INTERSTITIAL: InterstitialDefectConfiguration,
    }

    @classmethod
    def get_constructor(cls, defect_type: PointDefectTypeEnum) -> Type[PointDefectConfiguration]:
        try:
            return cls._type_to_configuration_map[defect_type]
        except KeyError:
            raise ValueError(f"Unsupported defect type: {defect_type}")


def create_defect_configuration(
    material: Union[Material, MaterialWithBuildMetadata],
    defect_type: PointDefectTypeEnum,
    coordinate: List[float],
    element: Optional[str] = None,
    placement_method: Union[
        VacancyPlacementMethodEnum,
        SubstitutionPlacementMethodEnum,
        InterstitialPlacementMethodEnum,
        AtomPlacementMethodEnum,
    ] = AtomPlacementMethodEnum.EXACT_COORDINATE,
    use_cartesian_coordinates: bool = False,
) -> PointDefectConfiguration:
    resolved_coordinate = resolve_coordinate(material, coordinate, placement_method, use_cartesian_coordinates)
    config_class = PointDefectConfigurationFactory.get_constructor(defect_type=defect_type)

    return config_class.from_parameters(crystal=material, coordinate=resolved_coordinate, element=element)
