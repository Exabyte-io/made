from typing import Union

from mat3ra.made.material import Material
from .configuration import PassivationConfiguration
from .builders import (
    SurfacePassivationBuilder,
    CoordinationBasedPassivationBuilder,
    SurfacePassivationBuilderParameters,
    CoordinationBasedPassivationBuilderParameters,
    CoordinationAnalyzer,
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
) -> set:
    """
    Get the unique coordination numbers for the provided passivation configuration as a set type.
        Considers the coordination threshold and shadowing radius from the builder parameters if provided.

    Args:
        configuration (PassivationConfiguration): The configuration object.
        builder_parameters (CoordinationBasedPassivationBuilderParameters): The builder parameters.

    Returns:
        set: The unique coordination numbers.
    """
    coordination_analyzer = CoordinationAnalyzer(cutoff=cutoff)
    return coordination_analyzer.get_unique_coordination_numbers(configuration.slab)
