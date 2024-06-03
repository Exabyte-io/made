from pydantic import BaseModel

from .termination_pair import TerminationPair
from ..slab import Termination
from ..slab.configuration import SlabConfiguration


class InterfaceConfiguration(BaseModel):
    film_configuration: SlabConfiguration
    substrate_configuration: SlabConfiguration
    film_termination: Termination
    substrate_termination: Termination
    distance_z: float = 3.0
    vacuum: float = 5.0

    @property
    def termination_pair(self):
        return TerminationPair(self.film_termination, self.substrate_termination)
