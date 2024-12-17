from typing import Union

from mat3ra.made.material import Material

from .builders import SlabBasedNanoparticleBuilder, ASEBasedNanoparticleBuilder
from .configuration import (
    ASEBasedNanoparticleConfiguration,
    SlabBasedNanoparticleConfiguration,
    SphereSlabBasedNanoparticleConfiguration,
)


def create_nanoparticle(
    configuration: Union[
        SlabBasedNanoparticleBuilder, SphereSlabBasedNanoparticleConfiguration, ASEBasedNanoparticleConfiguration
    ],
    builder: Union[SlabBasedNanoparticleBuilder, ASEBasedNanoparticleBuilder] = SlabBasedNanoparticleBuilder(),
) -> Material:
    return builder.get_material(configuration)
