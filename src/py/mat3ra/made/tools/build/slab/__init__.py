from mat3ra.made.material import Material
from .builders import SlabBuilder
from .configuration import SlabConfiguration
from .termination import Termination


def create_slab(
    configuration: SlabConfiguration,
) -> Material:
    builder = SlabBuilder()

    return builder.get_material(configuration)


def create_slab_if_not(material: Material, default_slab_configuration: SlabConfiguration) -> Material:
    slab = material
    if not slab.metadata or slab.metadata["build"]["configuration"]["type"] != SlabConfiguration.__name__:
        print("The material is not a slab. Creating a new slab...")
        slab = create_slab(default_slab_configuration)
    return slab
