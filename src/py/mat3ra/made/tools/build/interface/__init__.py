from typing import Union, List

from mat3ra.made.material import Material
from .builders import (
    SimpleInterfaceBuilder,
    SimpleInterfaceBuilderParameters,
    ZSLStrainMatchingParameters,
    ZSLStrainMatchingInterfaceBuilder,
    ZSLStrainMatchingInterfaceBuilderParameters,
)
from .configuration import InterfaceConfiguration


def create_interfaces(
    builder: Union[SimpleInterfaceBuilder, ZSLStrainMatchingInterfaceBuilder], configuration: InterfaceConfiguration
) -> List[Material]:
    return builder.get_materials(configuration)


def create_interface(
    builder: Union[SimpleInterfaceBuilder, ZSLStrainMatchingInterfaceBuilder],
    configuration: InterfaceConfiguration,
) -> Material:
    return builder.get_material(configuration)
