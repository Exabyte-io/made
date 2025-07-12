from typing import List

from ase.geometry import distance

from mat3ra.made.tools.build.defect.slab.helpers import recreate_slab_with_fractional_layers

from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.slab import SlabMaterialAnalyzer
from mat3ra.made.tools.analyze.crystal_site import CrystalSiteAnalyzer
from mat3ra.made.utils import get_atomic_coordinates_extremum
from mat3ra.made.tools.build.defect.enums import AdatomPlacementMethodEnum
from mat3ra.made.tools.build.slab.configurations import SlabStrainedSupercellConfiguration


# class SlabStackAnalyzer(SlabMaterialAnalyzer):
#     @property
#     def coordinate_in_added_component(self) -> List[float]:
#         """
#         Get the coordinate in the added component based on the slab configuration.
#         This method is overridden to provide a specific implementation for slab stacks.
#         """
#         # The coordinate in the added component is the same as in the slab
#         xy_coordinate = self.coordinate[:2]
#         z_coordinate = distance_z
#         return z_coordinate


class AdatomMaterialAnalyzer(SlabMaterialAnalyzer):
    distance_z: float
    placement_method: AdatomPlacementMethodEnum
    coordinate: List[float]  # Add coordinate property

    @property
    def distance_z_crystal(self) -> float:
        return self.material.basis.cell.convert_point_to_crystal([0, 0, self.distance_z])[2]

    @property
    def added_component_one_layer(self) -> Material:
        return recreate_slab_with_fractional_layers(self.material, 1)

    @property
    def coordinate_in_added_component_from_crystal_site(self) -> List[float]:
        crystal_site_analyzer = CrystalSiteAnalyzer(material=self.added_component_one_layer, coordinate=self.coordinate)
        return crystal_site_analyzer.closest_site_coordinate

    @property
    def coordinate_in_added_component(self) -> List[float]:
        if self.placement_method == AdatomPlacementMethodEnum.NEW_CRYSTAL_SITE:
            return self.coordinate_in_added_component_from_crystal_site
        else:
            return self.coordinate
