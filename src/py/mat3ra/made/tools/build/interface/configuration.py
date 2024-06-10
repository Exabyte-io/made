from mat3ra.code.entity import InMemoryEntity
from pydantic import BaseModel

from .termination_pair import TerminationPair
from ..slab import Termination
from ..slab.configuration import SlabConfiguration


class InterfaceConfiguration(BaseModel, InMemoryEntity):
    film_configuration: SlabConfiguration
    substrate_configuration: SlabConfiguration
    film_termination: Termination
    substrate_termination: Termination
    distance_z: float = 3.0
    vacuum: float = 5.0

    @property
    def termination_pair(self):
        return TerminationPair(self.film_termination, self.substrate_termination)

    @property
    def _json(self):
        return {
            "type": "InterfaceConfiguration",
            "film_configuration": self.film_configuration.to_json(),
            "substrate_configuration": self.substrate_configuration.to_json(),
            "film_termination": str(self.film_termination),
            "substrate_termination": str(self.substrate_termination),
            "distance_z": self.distance_z,
            "vacuum": self.vacuum,
        }
