from .termination_pair import TerminationPair
from ..slab.configuration import BaseSlabConfiguration, SlabConfiguration


class InterfaceConfiguration(BaseSlabConfiguration):
    film_configuration: SlabConfiguration
    substrate_configuration: SlabConfiguration
    film_termination: str
    substrate_termination: str
    distance_z: float = 3.0
    vacuum: float = 5.0

    @property
    def termination_pair(self):
        return TerminationPair((self.film_termination, self.substrate_termination))
