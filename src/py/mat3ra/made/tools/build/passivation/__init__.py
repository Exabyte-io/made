from typing import Union

from mat3ra.made.material import Material
from .configuration import PassivationConfiguration
from .builders import SurfacePassivationBuilder, UndercoordinationPassivationBuilder


def create_passivation(
    configuration: PassivationConfiguration,
    builder: Union[SurfacePassivationBuilder, UndercoordinationPassivationBuilder] = SurfacePassivationBuilder(),
) -> Material:
    return builder.get_material(configuration)
