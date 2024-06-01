from typing import List

from mat3ra.made.material import Material
from .builders import SlabBuilder, SlabSelectorParameters
from .configuration import SlabConfiguration


def get_terminations(configuration: SlabConfiguration) -> List[str]:
    return SlabBuilder().get_terminations(configuration)


def create_slab(configuration: SlabConfiguration, termination: str) -> Material:
    builder = SlabBuilder()
    return builder.get_material(configuration, selector_parameters=SlabSelectorParameters(termination=termination))
