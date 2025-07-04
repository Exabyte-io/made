from typing import Optional, List

from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.other import get_closest_site_id_from_coordinate
from mat3ra.made.tools.build.defect import PointDefectConfiguration
from mat3ra.made.tools.build.defect.enums import AtomPlacementMethodEnum, PointDefectTypeEnum
from mat3ra.made.tools.build.defect.point.configuration import (
    VacancyDefectConfiguration,
    SubstitutionalDefectConfiguration,
    InterstitialDefectConfiguration,
)


class PointDefectMaterialAnalyzer:
    """
    Analyzer for creating point defect configurations.

    Takes host material, defect type, coordinate, resolution method, and element
    to return appropriate point defect configurations.
    """

    def __init__(self):
        pass

    def create_point_defect_configuration(
        self,
        host_material: Material,
        defect_type: PointDefectTypeEnum,
        coordinate: Optional[List[float]] = None,
        resolution_method: AtomPlacementMethodEnum = AtomPlacementMethodEnum.COORDINATE,
        element: Optional[str] = None,
    ) -> PointDefectConfiguration:
        """
        Create a point defect configuration based on the provided parameters.

        Args:
            host_material: The host crystal material
            defect_type: Type of defect (vacancy, substitution, interstitial)
            coordinate: Optional coordinate for the defect
            resolution_method: Method to resolve coordinates
            element: Element for substitution/interstitial defects

        Returns:
            PointDefectConfiguration with appropriate crystal site configuration
        """

        if coordinate is None:
            coordinate = [0.0, 0.0, 0.0]

        resolved_coordinate = self._resolve_coordinate(
            host_material, coordinate, resolution_method, defect_type, element
        )

        if defect_type == PointDefectTypeEnum.VACANCY:
            return VacancyDefectConfiguration.from_parameters(
                host_material=host_material,
                coordinate=resolved_coordinate,
            )
        elif defect_type == PointDefectTypeEnum.SUBSTITUTION:
            if element is None:
                raise ValueError("Element must be specified for substitution defects")
            return SubstitutionalDefectConfiguration.from_parameters(
                host_material=host_material,
                coordinate=resolved_coordinate,
                element=element,
            )
        elif defect_type == PointDefectTypeEnum.INTERSTITIAL:
            if element is None:
                raise ValueError("Element must be specified for interstitial defects")
            return InterstitialDefectConfiguration.from_parameters(
                host_material=host_material,
                coordinate=resolved_coordinate,
                element=element,
            )
        else:
            raise ValueError(f"Unknown defect type: {defect_type}")

    def _resolve_coordinate(
        self,
        material: Material,
        coordinate: List[float],
        atom_placement: AtomPlacementMethodEnum,
        defect_type: PointDefectTypeEnum,
        element: Optional[str] = None,
    ) -> List[float]:

        if atom_placement == AtomPlacementMethodEnum.COORDINATE:
            return coordinate
        elif atom_placement == AtomPlacementMethodEnum.CLOSEST_SITE:
            return self._resolve_to_crystal_site(material, coordinate, element)
        elif atom_placement == AtomPlacementMethodEnum.NEW_CRYSTAL_SITE:
            return self._resolve_to_new_crystal_site(material, coordinate, element)
        elif atom_placement == AtomPlacementMethodEnum.EQUIDISTANT:
            return self._resolve_to_equidistant_site(material, coordinate)
        elif atom_placement == AtomPlacementMethodEnum.VORONOI_SITE:
            return self._resolve_to_voronoi_site(material, coordinate)
        else:
            raise ValueError(f"Unknown atom placement method: {atom_placement}")

    def _resolve_to_crystal_site(
        self, material: Material, coordinate: List[float], element: Optional[str] = None
    ) -> List[float]:
        site_id = get_closest_site_id_from_coordinate(material, coordinate)
        return material.coordinates_array[site_id]

    def _resolve_to_new_crystal_site(
        self, material: Material, coordinate: List[float], element: Optional[str] = None
    ) -> List[float]:
        # USe slab with additional layers to find crystal site in the next layer
        return coordinate

    def _resolve_to_equidistant_site(self, material: Material, coordinate: List[float]) -> List[float]:
        # For now, return the original coordinate
        # This could be enhanced with actual equidistant site calculation
        return coordinate

    def _resolve_to_voronoi_site(self, material: Material, coordinate: List[float]) -> List[float]:
        # For now, return the original coordinate
        # This could be enhanced with actual Voronoi site calculation
        return coordinate
