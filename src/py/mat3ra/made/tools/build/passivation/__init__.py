from typing import Union, List

from mat3ra.made.material import Material
from ...analyze.coordination import CoordinationAnalyzer
from .configuration import PassivationConfiguration
from .builders import (
    SurfacePassivationBuilder,
    CoordinationBasedPassivationBuilder,
    SurfacePassivationBuilderParameters,
    CoordinationBasedPassivationBuilderParameters,
)


def create_passivation(
    configuration: PassivationConfiguration,
    builder: Union[SurfacePassivationBuilder, CoordinationBasedPassivationBuilder, None] = None,
) -> Material:
    if builder is None:
        builder = SurfacePassivationBuilder(build_parameters=SurfacePassivationBuilderParameters())
    return builder.get_material(configuration)


def get_unique_coordination_numbers(
    configuration: PassivationConfiguration,
    cutoff: float = 3.0,
) -> List[int]:
    """
    Get the unique coordination numbers for the provided passivation configuration and cutoff radius.

    Args:
        configuration (PassivationConfiguration): The configuration object.
        cutoff (float): The cutoff radius for defining neighbors.
    Returns:
        set: The unique coordination numbers.
    """
    coordination_analyzer = CoordinationAnalyzer(cutoff=cutoff)
    return coordination_analyzer.get_unique_coordination_numbers(configuration.slab)
