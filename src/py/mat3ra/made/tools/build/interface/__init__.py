from typing import Union, List, Optional

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
    configuration: InterfaceConfiguration,
    builder: Optional[Union[SimpleInterfaceBuilder, ZSLStrainMatchingInterfaceBuilder]] = None,
) -> Material:
    if builder is None:
        builder = SimpleInterfaceBuilder(build_parameters=SimpleInterfaceBuilderParameters())
    return builder.get_material(configuration)
