from typing import List

import numpy as np

from mat3ra.made.material import Material

from ..interface.configuration import TwistedInterfaceConfiguration
from ..interface.builders import (
    CommensurateLatticeTwistedInterfaceBuilder,
    CommensurateLatticeTwistedInterfaceBuilderParameters,
)
from ..supercell import create_supercell
from ..utils import merge_two_materials_laterally


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
            phase_1_material_initial = create_supercell(item.configuration.film, matrix1.tolist())
            phase_2_material_initial = create_supercell(item.configuration.film, matrix2.tolist())

            interface = merge_two_materials_laterally(
                phase_1_material_initial,
                phase_2_material_initial,
                gap=item.configuration.gap,
                distance_tolerance=post_process_parameters.distance_tolerance,
            )
            grain_boundaries.append(interface)

        return grain_boundaries

    def _update_material_name(self, material: Material, configuration: SurfaceGrainBoundaryConfiguration) -> Material:
        material.name = f"Surface Grain Boundary ({configuration.twist_angle:.2f}°)"
        return material
