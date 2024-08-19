from mat3ra.made.material import Material

from .builders import NanoribbonBuilder
from .configuration import NanoribbonConfiguration


def build_nanoribbon(configuration: NanoribbonConfiguration) -> Material:
    builder = NanoribbonBuilder()
    return builder.get_material(configuration)
