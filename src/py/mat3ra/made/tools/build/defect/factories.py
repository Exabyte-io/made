from typing import Callable, List, Dict

from mat3ra.made.material import Material
from mat3ra.made.tools.build.defect.enums import PointDefectTypeEnum
from mat3ra.made.tools.build.defect.point.configuration import (
    VacancyDefectConfiguration,
    SubstitutionalDefectConfiguration,
    InterstitialDefectConfiguration,
    PointDefectConfiguration,
)

DEFECT_CONFIG_FACTORIES: Dict[PointDefectTypeEnum, Callable] = {
    PointDefectTypeEnum.VACANCY: lambda material, coordinate, element=None: VacancyDefectConfiguration.from_parameters(
        crystal=material, coordinate=coordinate
    ),
    PointDefectTypeEnum.SUBSTITUTION: lambda material, coordinate, element: SubstitutionalDefectConfiguration.from_parameters(
        crystal=material, coordinate=coordinate, element=element
    ),
    PointDefectTypeEnum.INTERSTITIAL: lambda material, coordinate, element: InterstitialDefectConfiguration.from_parameters(
        crystal=material, coordinate=coordinate, element=element
    ),
}


def create_defect_configuration(
    material: Material, defect_type: PointDefectTypeEnum, coordinate: List[float], element: str = None
) -> PointDefectConfiguration:
    if defect_type not in DEFECT_CONFIG_FACTORIES:
        raise ValueError(f"Unknown defect type: {defect_type}")

    factory = DEFECT_CONFIG_FACTORIES[defect_type]
    return factory(material, coordinate, element)
