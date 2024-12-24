from typing import List, Optional

from mat3ra.made.material import Material
from .builders import SlabBuilder, SlabSelectorParameters, PymatgenSlabGeneratorParameters
from .configuration import SlabConfiguration
from .termination import Termination


def get_terminations(
    configuration: SlabConfiguration, build_parameters: Optional[PymatgenSlabGeneratorParameters] = None
) -> List[Termination]:
    if build_parameters:
        return SlabBuilder(build_parameters=build_parameters).get_terminations(configuration)
    return SlabBuilder().get_terminations(configuration)


def create_slab(
    configuration: SlabConfiguration,
    termination: Optional[Termination] = None,
    build_parameters: Optional[PymatgenSlabGeneratorParameters] = None,
) -> Material:
    if build_parameters:
        builder = SlabBuilder(build_parameters=build_parameters)
    else:
        builder = SlabBuilder()
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
