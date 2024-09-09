from typing import Union

from mat3ra.made.material import Material
from .configuration import PassivationConfiguration
from .builders import (
    SurfacePassivationBuilder,
    UndercoordinationPassivationBuilder,
    SurfacePassivationBuilderParameters,
)


def create_passivation(
    configuration: PassivationConfiguration,
    builder: Union[SurfacePassivationBuilder, UndercoordinationPassivationBuilder, None] = None,
) -> Material:
    if builder is None:
        builder = SurfacePassivationBuilder(build_parameters=SurfacePassivationBuilderParameters())
    return builder.get_material(configuration)
