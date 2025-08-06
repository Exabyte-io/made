from typing import List

from mat3ra.made.material import Material
from mat3ra.made.tools.build.pristine_structures.two_dimensional.slab.builder import SlabBuilder

from ...build_components.entities.reusable.two_dimensional.slab_stack.helpers import (
    recreate_slab_with_fractional_layers,
)
from ...build_components.metadata import MaterialWithBuildMetadata
from .adatom_material_analyzer import AdatomMaterialAnalyzer
from .crystal_site_analyzer import CrystalSiteAnalyzer


class AdatomCrystalSiteMaterialAnalyzer(AdatomMaterialAnalyzer):
    DEFAULT_NUMBER_OF_LAYERS: float = 1

    @property
    def added_component_prototype(self) -> Material:
        # Recreate the slab with a single layer to ensure the adatom is placed correctly
        return recreate_slab_with_fractional_layers(self.material, self.DEFAULT_NUMBER_OF_LAYERS)

    @property
    def coordinate_in_added_component(self) -> List[float]:
        approximate_coordinate_3d = super().coordinate_in_added_component
        crystal_site_analyzer = CrystalSiteAnalyzer(
            material=self.added_component_prototype,
            coordinate=approximate_coordinate_3d,
        )
        return crystal_site_analyzer.closest_site_coordinate

    @property
    def slab_material_or_configuration_for_stacking(self) -> MaterialWithBuildMetadata:
        config = self.slab_configuration_with_no_vacuum
        params = self.build_parameters
        return SlabBuilder(build_parameters=params).get_material(config)
