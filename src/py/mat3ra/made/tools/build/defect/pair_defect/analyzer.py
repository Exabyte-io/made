from typing import List

from mat3ra.made.tools.analyze.crystal_site import CrystalSiteAnalyzer


class PairDefectMaterialAnalyzer(CrystalSiteAnalyzer):
    @property
    def primary_coordinate(self) -> List[float]:
        self.coordinate = self.primary_defect_configuration.merge_components[1].coordinate
        return super().closest_site_coordinate

    @property
    def secondary_coordinate(self) -> List[float]:
        self.coordinate = self.secondary_defect_configuration.merge_components[1].coordinate
        return super().closest_site_coordinate
