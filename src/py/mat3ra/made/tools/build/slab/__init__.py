from typing import List, Optional

from mat3ra.made.material import Material
from .builders import SlabBuilder, SlabSelectorParameters, PymatgenSlabGeneratorParameters
from .configuration import SlabConfiguration
from .termination import Termination

CACHED_BUILDER = None


def get_terminations(
    configuration: SlabConfiguration, build_parameters: Optional[PymatgenSlabGeneratorParameters] = None
) -> List[Termination]:
    global CACHED_BUILDER
    CACHED_BUILDER = SlabBuilder(build_parameters=build_parameters)
    return CACHED_BUILDER.get_terminations(configuration)


def create_slab(
    configuration: SlabConfiguration,
    termination: Optional[Termination] = None,
    build_parameters: Optional[PymatgenSlabGeneratorParameters] = None,
    use_cached_builder: bool = True,
) -> Material:
    builder = (
        CACHED_BUILDER if use_cached_builder and CACHED_BUILDER else SlabBuilder(build_parameters=build_parameters)
    )
    termination = termination or builder.get_terminations(configuration)[0]
    return builder.get_material(configuration, selector_parameters=SlabSelectorParameters(termination=termination))


def create_slab_if_not(material: Material, default_slab_configuration: SlabConfiguration) -> Material:
    """
    Create a slab from the material if it is not a slab already. Otherwise, return the material.

    Args:
        material (Material): The material to be checked.
        default_slab_configuration (SlabConfiguration): The default configuration to be used for creating a new slab.

    Returns:
        Material: The slab.
    """
    slab = material
    if not slab.metadata or slab.metadata["build"]["configuration"]["type"] != SlabConfiguration.__name__:
        print("The material is not a slab. Creating a new slab...")
        slab = create_slab(default_slab_configuration)
    return slab
