from .strained_configurations import SlabStrainedSupercellWithGapConfiguration
from ...vacuum.configuration import VacuumConfiguration


class SlabWithAdditionalLayersConfiguration(SlabStrainedSupercellWithGapConfiguration):
    type: str = "SlabWithAdditionalLayersConfiguration"
    number_of_additional_layers: float = 1.0

    @property
    def atomic_layers(self):
        return self.stack_components[0]

    @property
    def vacuum_configuration(self) -> VacuumConfiguration:
        return self.stack_components[1]
