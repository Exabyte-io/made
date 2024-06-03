from typing import List

from mat3ra.made.material import Material
from .builders import SlabBuilder, SlabSelectorParameters
from .configuration import SlabConfiguration
from .termination import Termination


def get_terminations(configuration: SlabConfiguration) -> List[Termination]:
    return SlabBuilder().get_terminations(configuration)


def create_slab(configuration: SlabConfiguration, termination: Termination) -> Material:
    builder = SlabBuilder()
    return builder.get_material(configuration, selector_parameters=SlabSelectorParameters(termination=termination))
