from typing import List

from mat3ra.esse.models.materials_category_components.entities.core.zero_dimensional.atom import AtomSchema
from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.crystal_site import CrystalSiteAnalyzer
from mat3ra.made.tools.analyze.slab import SlabMaterialAnalyzer
from mat3ra.made.tools.build import MaterialWithBuildMetadata
from mat3ra.made.tools.build.defect.point.builders import AtomAtCoordinateBuilder, AtomAtCoordinateConfiguration
from mat3ra.made.tools.build.defect.slab.helpers import recreate_slab_with_fractional_layers
from mat3ra.made.tools.build.vacuum.builders import VacuumBuilder
from mat3ra.made.tools.build.vacuum.configuration import VacuumConfiguration


class AdatomMaterialAnalyzer(SlabMaterialAnalyzer):
    distance_z: float
    coordinate_2d: List[float]  # Add coordinate property
    element: str

    @property
    def added_component_height(self) -> float:
        return self.distance_z * 2

    @property
    def added_component_prototype(self) -> MaterialWithBuildMetadata:
        vacuum_configuration = VacuumConfiguration(crystal=self.material, size=self.added_component_height)
        vacuum_material = VacuumBuilder().get_material(vacuum_configuration)
        return vacuum_material

    @property
    def coordinate_in_added_component(self) -> List[float]:
        coordinate_3d = self.coordinate_2d + [0]
        coordinate_3d_cartesian = self.material.basis.cell.convert_point_to_cartesian(coordinate_3d)
        coordinate_3d_cartesian[2] = self.distance_z
        coordinate_3d_crystal = self.added_component_prototype.basis.cell.convert_point_to_crystal(
            coordinate_3d_cartesian
        )
        return coordinate_3d_crystal

    @property
    def added_component(self) -> MaterialWithBuildMetadata:
        atom_configuration = AtomAtCoordinateConfiguration(
            crystal=self.added_component_prototype,
            element=AtomSchema(chemical_element=self.element),
            coordinate=self.coordinate_in_added_component,
        )
        return AtomAtCoordinateBuilder().get_material(atom_configuration)


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
