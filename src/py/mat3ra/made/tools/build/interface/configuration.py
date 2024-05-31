from src.py.mat3ra.made.material import Material
from .termination_pair import TerminationPair
from ..slab.configuration import BaseSlabConfiguration, SlabConfiguration
from .builders import SimpleInterfaceBuilder, SimpleInterfaceBuilderParameters


class InterfaceConfiguration(BaseSlabConfiguration):
    film_configuration: SlabConfiguration
    substrate_configuration: SlabConfiguration
    film_termination: str
    substrate_termination: str
    termination_pair: TerminationPair
    distance_z: float = 3.0
    vacuum: float = 5.0

    def __init__(
        self,
        film_configuration: SlabConfiguration,
        substrate_configuration: SlabConfiguration,
        film_termination: str,
        substrate_termination: str,
        distance_z: float = 3.0,
        vacuum: float = 5.0,
    ):
        super().__init__()
        self.film_configuration = film_configuration
        self.substrate_configuration = substrate_configuration
        self.film_termination = film_termination
        self.substrate_termination = substrate_termination
        self.termination_pair = TerminationPair((film_termination, substrate_termination))
        self.distance_z: float = distance_z
        self.vacuum: float = vacuum
        self.__builder = SimpleInterfaceBuilder(build_parameters=SimpleInterfaceBuilderParameters(scale_film=False))

    @property
    def bulk(self):
        # TODO: implement utils to remove vacuum
        return self.get_interface()

    def get_interface(self) -> Material:
        return self.__builder.get_material(self)
