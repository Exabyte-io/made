from typing import Optional, List

from .. import BaseConfiguration
from ..slab.configuration import SlabConfiguration
from ..slab.termination import Termination
from ..interface.configuration import TwistedInterfaceConfiguration


class SlabGrainBoundaryConfiguration(BaseConfiguration):
    """
    Configuration for a grain boundary between two phases with different surfaces facing each other.

    Attributes:
        phase_1_configuration (SlabConfiguration): The configuration of the first phase.
        phase_2_configuration (SlabConfiguration): The configuration of the second phase.
        phase_1_termination (Termination): The termination of the first phase.
        phase_2_termination (Termination): The termination of the second phase.
        gap (float): The gap between the two phases, in Angstroms.
        slab_configuration (SlabConfiguration): The configuration of the grain boundary slab.
        slab_termination (Optional[Termination]): The termination of the grain boundary slab.
    """

    phase_1_configuration: SlabConfiguration
    phase_2_configuration: SlabConfiguration
    phase_1_termination: Termination
    phase_2_termination: Termination
    gap: float = 3.0
    slab_configuration: SlabConfiguration
    slab_termination: Optional[Termination] = None

    @property
    def _json(self):
        return {
            "type": self.__class__.__name__,
            "phase_1_configuration": self.phase_1_configuration.to_json(),
            "phase_2_configuration": self.phase_2_configuration.to_json(),
            "phase_1_termination": str(self.phase_1_termination),
            "phase_2_termination": str(self.phase_2_termination),
            "gap": self.gap,
            "slab_configuration": self.slab_configuration.to_json(),
        }


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
