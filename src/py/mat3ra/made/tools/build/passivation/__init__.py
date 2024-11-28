from typing import Union, Optional

from mat3ra.made.material import Material
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
    builder_parameters: Optional[CoordinationBasedPassivationBuilderParameters] = None,
) -> set:
    """
    Get the unique coordination numbers for the provided configuration as a set type.

    Args:
        configuration (PassivationConfiguration): The configuration object.
        builder_parameters (CoordinationBasedPassivationBuilderParameters): The builder parameters.
    """
    if builder_parameters is None:
        builder_parameters = CoordinationBasedPassivationBuilderParameters()
    return CoordinationBasedPassivationBuilder(builder_parameters).get_unique_coordination_numbers(
        material=configuration.slab
    )
