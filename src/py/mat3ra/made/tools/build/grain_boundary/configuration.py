from typing import Optional

from .. import BaseConfiguration
from ..slab.configuration import SlabConfiguration
from ..slab.termination import Termination


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
