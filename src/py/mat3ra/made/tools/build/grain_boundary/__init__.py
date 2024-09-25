from typing import List

from mat3ra.made.tools.build import BaseConfiguration

from mat3ra.made.material import Material

from mat3ra.made.tools.build.interface.builders import CommensurateLatticeTwistedInterfaceBuilder
from mat3ra.made.tools.build.supercell import create_supercell
from mat3ra.made.tools.build.utils import merge_materials
from mat3ra.made.tools.modify import translate_by_vector, filter_by_box


class SurfaceGrainBoundaryConfiguration(BaseConfiguration):
    """
    Configuration for creating a surface grain boundary.

    Args:
        material (Material): The material to create two phases of the grain boundary.
        angle (float): The angle between the orientations of the two phases.
        gap (float): The gap between the two phases.
    """

    material: Material
    angle: float = 0.0
    gap: float = 0.0

    @property
    def _json(self):
        return {
            "type": self.get_cls_name(),
            "material": self.material.to_json(),
            "angle": self.angle,
        }


class SurfaceGrainBoundaryBuilder(CommensurateLatticeTwistedInterfaceBuilder):
    def _post_process(
        self,
        items: List[CommensurateLatticeTwistedInterfaceBuilder._GeneratedItemType],
        post_process_parameters=None,
    ) -> List[Material]:
        grain_boundary_materials = []
        for item in items:
            phase_1_material_initial = create_supercell(item.configuration.film, item.matrix1.tolist())
            phase_1_material_doubled = create_supercell(phase_1_material_initial, scaling_factors=[2, 1, 1])
            phase_1_material = filter_by_box(phase_1_material_doubled, [0, 0, 0], [0.5, 1, 1])
            phase_2_material_initial = create_supercell(item.configuration.substrate, item.matrix2.tolist())
            phase_2_material_doubled = create_supercell(phase_2_material_initial, scaling_factors=[2, 1, 1])
            phase_2_material = filter_by_box(phase_2_material_doubled, [0.5, 0, 0], [1, 1, 1])
            interface = merge_materials([phase_1_material, phase_2_material])
            interface.metadata["actual_twist_angle"] = item.angle
            grain_boundary_materials.append(interface)
        return grain_boundary_materials
