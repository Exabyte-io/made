from typing import Any, Generic, Optional, Type, TypeVar

from mat3ra.made.tools.build_components.metadata import MaterialWithBuildMetadata
from pydantic import BaseModel

from .base_configuration_pydantic import BaseConfigurationPydantic
from .build_parameters import BaseBuilderParameters

TypeConfiguration = TypeVar("TypeConfiguration", bound=BaseConfigurationPydantic)
TypeBuildParameters = TypeVar("TypeBuildParameters", bound=BaseBuilderParameters)


class BaseSingleBuilder(BaseModel, Generic[TypeConfiguration, TypeBuildParameters]):
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

    model_config = {"arbitrary_types_allowed": True}

    build_parameters: Optional[TypeBuildParameters] = None
    _BuildParametersType: Type[TypeBuildParameters] = BaseBuilderParameters
    _DefaultBuildParameters: TypeBuildParameters = BaseBuilderParameters()

    _ConfigurationType: Type[TypeConfiguration] = BaseConfigurationPydantic
    _PostProcessParametersType: Any = None

    def __init__(self, build_parameters: Optional[TypeBuildParameters] = None):
        super().__init__()
        self.build_parameters = build_parameters if build_parameters is not None else self._DefaultBuildParameters

    def _generate(self, configuration: TypeConfiguration) -> MaterialWithBuildMetadata:
        raise NotImplementedError

    def _post_process(
        self,
        item: MaterialWithBuildMetadata,
        post_process_parameters: Optional[Any],
        configuration: Optional[TypeConfiguration] = None,
    ) -> MaterialWithBuildMetadata:
        return item

    def _finalize(
        self, material: MaterialWithBuildMetadata, configuration: TypeConfiguration
    ) -> MaterialWithBuildMetadata:
        material_with_metadata = self._update_material_metadata(material, configuration)
        return self._update_material_name(material_with_metadata, configuration)

    def get_material(
        self,
        configuration: TypeConfiguration,
        post_process_parameters: Optional[Any] = None,
    ) -> MaterialWithBuildMetadata:
        generated_item = self._generate(configuration)
        material = self._post_process(generated_item, post_process_parameters, configuration)
        finalized_material = self._finalize(material, configuration)
        return finalized_material

    def _update_material_name(
        self, material: MaterialWithBuildMetadata, configuration: TypeConfiguration
    ) -> MaterialWithBuildMetadata:
        return material

    def _update_material_metadata(
        self, material: MaterialWithBuildMetadata, configuration: TypeConfiguration
    ) -> MaterialWithBuildMetadata:
        material.metadata.add_build_metadata_step(configuration=configuration, build_parameters=self.build_parameters)
        return material
