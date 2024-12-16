from mat3ra.made.material import Material

from .builders import NanoparticleBuilder
from .configuration import ASENanoparticleConfiguration


def create_nanoparticle(configuration: ASENanoparticleConfiguration) -> Material:
    builder = NanoparticleBuilder()
    return builder.get_material(configuration)
