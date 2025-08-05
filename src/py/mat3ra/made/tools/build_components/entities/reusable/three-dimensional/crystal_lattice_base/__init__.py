from typing import TypeVar

from .base_configuration_pydantic import BaseConfigurationPydantic
from .build_parameters import BaseBuilderParameters
from .base_single_builder import BaseSingleBuilder, TypeConfiguration, TypeBuildParameters

# Re-export the TypeVar for backwards compatibility
BaseConfigurationPydanticChild = TypeVar("BaseConfigurationPydanticChild", bound="BaseConfigurationPydantic")

__all__ = [
    "BaseConfigurationPydantic",
    "BaseBuilderParameters",
    "BaseSingleBuilder",
    "TypeConfiguration",
    "TypeBuildParameters",
    "BaseConfigurationPydanticChild",
]
