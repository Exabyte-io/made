from mat3ra.made.material import Material

from .builders import NanoparticleBuilder
from .configuration import NanoparticleConfiguration


def create_nanoparticle(configuration: NanoparticleConfiguration) -> Material:
    builder = NanoparticleBuilder()
    return builder.get_material(configuration)
