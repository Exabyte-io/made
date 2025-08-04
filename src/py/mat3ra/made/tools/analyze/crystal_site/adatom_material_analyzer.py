from typing import List

from mat3ra.esse.models.materials_category_components.entities.core.zero_dimensional.atom import AtomSchema
from mat3ra.made.tools.analyze.slab import SlabMaterialAnalyzer
from mat3ra.made.tools.build import MaterialWithBuildMetadata
from mat3ra.made.tools.build.defect.point.atom_at_coordinate.builder import AtomAtCoordinateBuilder
from mat3ra.made.tools.build.defect.point.atom_at_coordinate.configuration import AtomAtCoordinateConfiguration
from mat3ra.made.tools.build.vacuum.builders import VacuumBuilder
from mat3ra.made.tools.build.vacuum.configuration import VacuumConfiguration


class AdatomMaterialAnalyzer(SlabMaterialAnalyzer):
    distance_z: float
    coordinate_2d: List[float]  # Add coordinate property
    element: str

    @property
    def added_component_height(self) -> float:
        return self.layer_thickness

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

    @property
    def slab_material_or_configuration_for_stacking(self) -> MaterialWithBuildMetadata:
        return self._slab_with_no_gap
