from typing import List

from mat3ra.made.tools.build import BaseConfiguration

from mat3ra.made.material import Material

from mat3ra.made.tools.build.interface.builders import CommensurateLatticeTwistedInterfaceBuilder
from mat3ra.made.tools.build.supercell import create_supercell
from mat3ra.made.tools.build.utils import merge_materials
from mat3ra.made.tools.modify import filter_by_box, add_vacuum_sides, translate_by_vector

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

            new_lattice_vectors_1 = phase_1_material.lattice.vector_arrays
            new_lattice_vectors_1[0][0] += item.configuration.gap

            new_lattice_vectors_2 = phase_2_material.lattice.vector_arrays
            new_lattice_vectors_2[0][0] += item.configuration.gap
            phase_1_material.set_new_lattice_vectors(
                lattice_vector1=new_lattice_vectors_1[0],
                lattice_vector2=new_lattice_vectors_1[1],
                lattice_vector3=new_lattice_vectors_1[2],
            )
            phase_2_material.set_new_lattice_vectors(
                lattice_vector1=new_lattice_vectors_2[0],
                lattice_vector2=new_lattice_vectors_2[1],
                lattice_vector3=new_lattice_vectors_2[2],
            )

            phase_2_material = translate_by_vector(
                phase_2_material, [item.configuration.gap / 2, 0, 0], use_cartesian_coordinates=True
            )

            interface = merge_materials([phase_1_material, phase_2_material])
            grain_boundaries.append(interface)

        return grain_boundaries

    def _update_material_name(self, material: Material, configuration: SurfaceGrainBoundaryConfiguration) -> Material:
        material.name = f"Surface Grain Boundary ({configuration.twist_angle:.2f}°)"
        return material
