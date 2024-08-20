from mat3ra.made.material import Material

from .builders import NanoribbonBuilder
from .configuration import NanoribbonConfiguration


def create_nanoribbon(configuration: NanoribbonConfiguration) -> Material:
    builder = NanoribbonBuilder()
    return builder.get_material(configuration)
