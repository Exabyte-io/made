from typing import Optional

from .. import BaseConfiguration
from ..slab.configuration import SlabConfiguration
from ..slab.termination import Termination


class GrainBoundaryConfiguration(BaseConfiguration):
    """
    Configuration for a grain boundary in a slab material.

    Attributes:
        phase_1_configuration: SlabConfiguration
        phase_2_configuration: SlabConfiguration
        phase_1_termination: Termination
        phase_2_termination: Termination
        gap: float
        slab_configuration: SlabConfiguration
        slab_termination: Termination

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
            "type": "GrainBoundaryConfiguration",
            "phase_1_configuration": self.phase_1_configuration.to_json(),
            "phase_2_configuration": self.phase_2_configuration.to_json(),
            "phase_1_termination": str(self.phase_1_termination),
            "phase_2_termination": str(self.phase_2_termination),
            "gap": self.gap,
            "slab_configuration": self.slab_configuration.to_json(),
        }
