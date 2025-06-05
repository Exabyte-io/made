from typing import List, Optional

from mat3ra.made.material import Material
from .builders import SlabBuilder, SlabBuilderParameters
from .configuration import SlabConfiguration
from .termination import Termination

CACHED_BUILDER = None

#
# def get_terminations(
#     configuration: SlabConfiguration, build_parameters: Optional[SlabBuilderParameters] = None
# ) -> List[Termination]:
#     global CACHED_BUILDER
#     CACHED_BUILDER = SlabBuilder(build_parameters=build_parameters)
#     return CACHED_BUILDER.get_terminations(configuration)


def create_slab(
    configuration: SlabConfiguration,
    build_parameters: Optional[SlabBuilderParameters] = None,
    use_cached_builder: bool = False,
) -> Material:
    builder = (
        CACHED_BUILDER if use_cached_builder and CACHED_BUILDER else SlabBuilder(build_parameters=build_parameters)
    )
    return builder.get_material(configuration)


def create_slab_if_not(material: Material, default_slab_configuration: SlabConfiguration) -> Material:
    slab = material
    if not slab.metadata or slab.metadata["build"]["configuration"]["type"] != SlabConfiguration.__name__:
        print("The material is not a slab. Creating a new slab...")
        slab = create_slab(default_slab_configuration)
    return slab
