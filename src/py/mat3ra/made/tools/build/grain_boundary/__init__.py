from typing import List

from mat3ra.made.tools.build import BaseConfiguration

from mat3ra.made.material import Material

from mat3ra.made.tools.build.interface.builders import CommensurateLatticeTwistedInterfaceBuilder
from mat3ra.made.tools.build.supercell import create_supercell
from mat3ra.made.tools.build.utils import merge_materials
from mat3ra.made.tools.modify import filter_by_box


from mat3ra.made.tools.build.interface.configuration import TwistedInterfaceConfiguration


class SurfaceGrainBoundaryConfiguration(TwistedInterfaceConfiguration):
    """
    Configuration for creating a surface grain boundary.

    Args:
        gap (float): The gap between the two phases.
    """

    gap: float = 0.0

    @property
    def _json(self):
        return {
            "type": self.get_cls_name(),
            "gap": self.gap,
        }


class SurfaceGrainBoundaryBuilder(CommensurateLatticeTwistedInterfaceBuilder):
    _ConfigurationType: type(SurfaceGrainBoundaryConfiguration) = SurfaceGrainBoundaryConfiguration  # type: ignore

    def _post_process(self, items: List[Material], post_process_parameters=None) -> List[Material]:
        grain_boundaries = []
        for item in items:
            phase_1_material_initial = create_supercell(item.configuration.film, item.matrix1.tolist())
            phase_1_material_doubled = create_supercell(phase_1_material_initial, scaling_factor=[2, 1, 1])
            phase_1_material = filter_by_box(phase_1_material_doubled, [0, 0, 0], [0.5, 1, 1])
            phase_2_material_initial = create_supercell(item.configuration.film, item.matrix2.tolist())
            phase_2_material_doubled = create_supercell(phase_2_material_initial, scaling_factor=[2, 1, 1])
            phase_2_material = filter_by_box(phase_2_material_doubled, [0.5, 0, 0], [1, 1, 1])

            interface = merge_materials([phase_1_material, phase_2_material])
            grain_boundaries.append(interface)

        return grain_boundaries

    def _update_material_name(self, material: Material, configuration: SurfaceGrainBoundaryConfiguration) -> Material:
        material.name = f"Surface Grain Boundary ({configuration.twist_angle:.2f}°)"
        return material
