from typing import List, Optional, Any, TypeVar

from mat3ra.code.entity import InMemoryEntityPydantic, InMemoryEntity
from pydantic import BaseModel

from .metadata import MaterialBuildMetadata, BuildMetadata, MaterialWithBuildMetadata
from ...material import Material

BaseConfigurationPydanticChild = TypeVar("BaseConfigurationPydanticChild", bound="BaseConfigurationPydantic")


class BaseConfiguration(BaseModel, InMemoryEntity):
    """
    Base class for material build configurations.
    This class provides an interface for defining the configuration parameters.

    The class is designed to be subclassed and the subclass should define the following attributes:
    - `_json`: The JSON representation of the configuration.
    """

    class Config:
        arbitrary_types_allowed = True

    @property
    def _json(self):
        raise NotImplementedError


class BaseConfigurationPydantic(InMemoryEntityPydantic):
    pass


class BaseSelectorParameters(BaseModel):
    default_index: int = 0


class BaseBuilderParameters(InMemoryEntityPydantic):
    pass


class BaseSingleBuilder(BaseModel):
    """
    Base class for material builders.
    This class provides an interface for generating materials and getter functions.
    The builder is meant as a description of the process, while its functions require a
    "Configuration" class instance to perform the generation.

    The class is designed to be subclassed and the subclass should implement the following methods:

    - `_generate`: Generate the material items, possibly using third-party tools/implementation for items.
    - `_sort`: Sort the items.
    - `_select`: Select a subset of the items.
    - `_post_process`: Post-process the items to convert them to materials (Material class).
    - `_finalize`: Finalize the materials.

    The subclass should also define the following attributes:

    - `_BuildParametersType`: The data structure model for the build parameters.
    - `_DefaultBuildParameters`: The default build parameters.
    - `_ConfigurationType`: The data structure model for the Configuration used during the build.
    - `_GeneratedItemType`: The type of the generated item.
    - `_SelectorParametersType`: The data structure model for the selector parameters.
    - `_PostProcessParametersType`: The data structure model for the post-process parameters.
    """

    class Config:
        arbitrary_types_allowed = True

    build_parameters: Any = None
    _BuildParametersType: Any = None
    _DefaultBuildParameters: Any = None

    _ConfigurationType: Any = Any
    _PostProcessParametersType: Any = None

    def __init__(self, build_parameters: _BuildParametersType = None):
        super().__init__(build_parameters=build_parameters)
        self.build_parameters = build_parameters if build_parameters is not None else self._DefaultBuildParameters

    def _generate(self, configuration: _ConfigurationType) -> MaterialWithBuildMetadata:
        raise NotImplementedError

    def _post_process(
        self, item: MaterialWithBuildMetadata, post_process_parameters: Optional[_PostProcessParametersType]
    ) -> MaterialWithBuildMetadata:
        return item

    def _finalize(
        self, material: MaterialWithBuildMetadata, configuration: _ConfigurationType
    ) -> MaterialWithBuildMetadata:
        material_with_metadata = self._update_material_metadata(material, configuration)
        return self._update_material_name(material_with_metadata, configuration)

    def get_material(
        self,
        configuration: _ConfigurationType,
        post_process_parameters: Optional[_PostProcessParametersType] = None,
    ) -> MaterialWithBuildMetadata:
        generated_item = self._generate(configuration)
        material = self._post_process(generated_item, post_process_parameters)
        finalized_material = self._finalize(material, configuration)
        return finalized_material

    def _update_material_name(self, material, configuration) -> MaterialWithBuildMetadata:
        return material

    def _update_material_metadata(
        self, material: MaterialWithBuildMetadata, configuration
    ) -> MaterialWithBuildMetadata:
        material.metadata.add_build_metadata_step(configuration=configuration, build_parameters=self.build_parameters)
        return material
