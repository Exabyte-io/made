from typing import Union

from mat3ra.made.material import Material

from .builders import NanoparticleBuilder, ASEBasedNanoparticleBuilder
from .configuration import ASENanoparticleConfiguration, NanoparticleConfiguration


def create_nanoparticle(
    configuration: Union[NanoparticleConfiguration, ASENanoparticleConfiguration],
    builder: Union[NanoparticleBuilder, ASEBasedNanoparticleBuilder] = NanoparticleBuilder(),
) -> Material:
    return builder.get_material(configuration)
