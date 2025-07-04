from typing import List, Optional

from mat3ra.made.tools.analyze import BaseMaterialAnalyzer
from mat3ra.made.tools.analyze.other import get_closest_site_id_from_coordinate
from mat3ra.made.tools.build.defect.enums import AtomPlacementMethodEnum, PointDefectTypeEnum
from mat3ra.made.tools.build.defect.point.configuration import (
    InterstitialDefectConfiguration,
    PointDefectConfiguration,
    SubstitutionalDefectConfiguration,
    VacancyDefectConfiguration,
)


class PointDefectAnalyzer(BaseMaterialAnalyzer):
    resolution_method: AtomPlacementMethodEnum = AtomPlacementMethodEnum.COORDINATE

    def _resolve_coordinate(self, coordinate: List[float], element: Optional[str] = None) -> List[float]:
        if self.resolution_method == AtomPlacementMethodEnum.COORDINATE:
            return coordinate
        elif self.resolution_method == AtomPlacementMethodEnum.CLOSEST_SITE:
            return self._resolve_to_crystal_site(coordinate, element)
        elif self.resolution_method == AtomPlacementMethodEnum.NEW_CRYSTAL_SITE:
            return self._resolve_to_new_crystal_site(coordinate, element)
        elif self.resolution_method == AtomPlacementMethodEnum.EQUIDISTANT:
            return self._resolve_to_equidistant_site(coordinate)
        elif self.resolution_method == AtomPlacementMethodEnum.VORONOI_SITE:
            return self._resolve_to_voronoi_site(coordinate)
        else:
            raise ValueError(f"Unknown atom placement method: {self.resolution_method}")

    def _resolve_to_crystal_site(self, coordinate: List[float], element: Optional[str] = None) -> List[float]:
        site_id = get_closest_site_id_from_coordinate(self.material, coordinate)
        return self.material.coordinates_array[site_id]

    def _resolve_to_new_crystal_site(self, coordinate: List[float], element: Optional[str] = None) -> List[float]:
        return coordinate

    def _resolve_to_equidistant_site(self, coordinate: List[float]) -> List[float]:
        return coordinate

    def _resolve_to_voronoi_site(self, coordinate: List[float]) -> List[float]:
        return coordinate

    def get_configuration(
        self,
        defect_type: PointDefectTypeEnum,
        coordinate: Optional[List[float]] = None,
        element: Optional[str] = None,
    ) -> PointDefectConfiguration:
        if coordinate is None:
            coordinate = [0.0, 0.0, 0.0]

        resolved_coordinate = self._resolve_coordinate(coordinate, element)

        if defect_type == PointDefectTypeEnum.VACANCY:
            return VacancyDefectConfiguration.from_parameters(crystal=self.material, coordinate=resolved_coordinate)
        elif defect_type == PointDefectTypeEnum.SUBSTITUTION:
            if element is None:
                raise ValueError("Element must be specified for substitution defects")
            return SubstitutionalDefectConfiguration.from_parameters(
                crystal=self.material, coordinate=resolved_coordinate, element=element
            )
        elif defect_type == PointDefectTypeEnum.INTERSTITIAL:
            if element is None:
                raise ValueError("Element must be specified for interstitial defects")
            return InterstitialDefectConfiguration.from_parameters(
                crystal=self.material, coordinate=resolved_coordinate, element=element
            )
        else:
            raise ValueError(f"Unknown defect type: {defect_type}")
