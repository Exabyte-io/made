from typing import Union

from mat3ra.made.material import Material

from .builders import NanoparticleBuilder, ASEBasedNanoparticleBuilder
from .configuration import ASEBasedNanoparticleConfiguration, NanoparticleConfiguration


def create_nanoparticle(
    configuration: Union[NanoparticleConfiguration, ASEBasedNanoparticleConfiguration],
    builder: Union[NanoparticleBuilder, ASEBasedNanoparticleBuilder] = NanoparticleBuilder(),
) -> Material:
    return builder.get_material(configuration)
