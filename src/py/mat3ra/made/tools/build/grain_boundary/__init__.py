from typing import List

import numpy as np
from mat3ra.made.tools.build import BaseConfiguration

from mat3ra.made.material import Material

from ..interface.builders import (
    CommensurateLatticeTwistedInterfaceBuilder,
    CommensurateLatticeTwistedInterfaceBuilderParameters,
)
from ..supercell import create_supercell
from ..utils import merge_materials
from ...modify import filter_by_box, add_vacuum_sides, translate_by_vector

from ..interface.configuration import TwistedInterfaceConfiguration


class SurfaceGrainBoundaryConfiguration(TwistedInterfaceConfiguration):
    """
    Configuration for creating a surface grain boundary.

    Args:
        gap (float): The gap between the two phases.
    """

    gap: float = 0.0
    xy_supercell_matrix: List[List[int]] = [[1, 0], [0, 1]]

    @property
    def _json(self):
        return {
            "type": self.get_cls_name(),
            "gap": self.gap,
            "xy_supercell_matrix": self.xy_supercell_matrix,
        }


class SurfaceGrainBoundaryBuilderParameters(CommensurateLatticeTwistedInterfaceBuilderParameters):
    """
    Parameters for creating a surface grain boundary.

    Args:
        distance_tolerance (float): The distance tolerance to remove atoms that are too close, in angstroms.
    """

    distance_tolerance: float = 0.1


class SurfaceGrainBoundaryBuilder(CommensurateLatticeTwistedInterfaceBuilder):
    _ConfigurationType: type(SurfaceGrainBoundaryConfiguration) = SurfaceGrainBoundaryConfiguration  # type: ignore
    _BuildParametersType = SurfaceGrainBoundaryBuilderParameters
    _DefaultBuildParameters = SurfaceGrainBoundaryBuilderParameters()

    def _post_process(self, items: List[Material], post_process_parameters=None) -> List[Material]:
        grain_boundaries = []
        for item in items:
            matrix1 = np.dot(np.array(item.configuration.xy_supercell_matrix), item.matrix1)
            matrix2 = np.dot(np.array(item.configuration.xy_supercell_matrix), item.matrix2)
            # TODO: isolate placing side by side with a gap into a separate function (in interface.utils folder?)
            phase_1_material_initial = create_supercell(item.configuration.film, matrix1.tolist())
            phase_1_material_doubled = create_supercell(phase_1_material_initial, scaling_factor=[2, 1, 1])
            phase_1_material = filter_by_box(phase_1_material_doubled, [0, 0, 0], [0.5, 1, 1])
            phase_2_material_initial = create_supercell(item.configuration.film, matrix2.tolist())
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

            interface = merge_materials(
                [phase_1_material, phase_2_material], distance_tolerance=self.build_parameters.distance_tolerance
            )
            grain_boundaries.append(interface)

        return grain_boundaries

    def _update_material_name(self, material: Material, configuration: SurfaceGrainBoundaryConfiguration) -> Material:
        material.name = f"Surface Grain Boundary ({configuration.twist_angle:.2f}°)"
        return material
