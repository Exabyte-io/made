from typing import Union

from mat3ra.made.material import Material
from .configuration import PassivationConfiguration
from .builders import (
    SurfacePassivationBuilder,
    CoordinationBasedPassivationBuilder,
    SurfacePassivationBuilderParameters,
)


def create_passivation(
    configuration: PassivationConfiguration,
    builder: Union[SurfacePassivationBuilder, CoordinationBasedPassivationBuilder, None] = None,
) -> Material:
    if builder is None:
        builder = SurfacePassivationBuilder(build_parameters=SurfacePassivationBuilderParameters())
    return builder.get_material(configuration)


def get_unique_coordination_numbers(configuration: PassivationConfiguration) -> set:
    """
    Get the unique coordination numbers for the provided configuration as a set type.
    """
    return CoordinationBasedPassivationBuilder().get_unique_coordination_numbers(material=configuration.material)
