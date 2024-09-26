from typing import List

from ..interface.configuration import TwistedInterfaceConfiguration


class SurfaceGrainBoundaryConfiguration(TwistedInterfaceConfiguration):
    """
    Configuration for creating a surface grain boundary.

    Args:
        gap (float): The gap between the two phases.
        xy_supercell_matrix (List[List[int]]): The supercell matrix to apply for both phases.
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
