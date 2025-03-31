from typing import List, Optional, Any

from mat3ra.code.entity import InMemoryEntity
from pydantic import BaseModel

from ...material import Material


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


class BaseSelectorParameters(BaseModel):
    default_index: int = 0


class BaseBuilderParameters(BaseModel):
    class Config:
        arbitrary_types_allowed = True


class BaseBuilder(BaseModel):
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
    _GeneratedItemType: Any = Any
    _SelectorParametersType: Any = None
    _PostProcessParametersType: Any = None
    selector_parameters: Any = BaseSelectorParameters()

    def __init__(self, build_parameters: _BuildParametersType = None):
        super().__init__(build_parameters=build_parameters)
        self.build_parameters = build_parameters or self._DefaultBuildParameters
        self.__generated_items: List[List[BaseBuilder._GeneratedItemType]] = []
        self.__configurations: List[BaseBuilder._ConfigurationType] = []

    def _generate(self, configuration: _ConfigurationType) -> List[_GeneratedItemType]:
        return []

    def _generate_or_get_from_cache(self, configuration: _ConfigurationType) -> List[_GeneratedItemType]:
        if configuration not in self.__configurations:
            self.__configurations.append(configuration)
            self.__generated_items.append(self._generate(configuration))
        return self.__generated_items[self.__configurations.index(configuration)]

    def _sort(self, items: List[_GeneratedItemType]) -> List[_GeneratedItemType]:
        return items

    def _select(
        self, items: List[_GeneratedItemType], selector_parameters: Optional[_SelectorParametersType]
    ) -> List[_GeneratedItemType]:
        return items

    def _post_process(
        self, items: List[_GeneratedItemType], post_process_parameters: Optional[_PostProcessParametersType]
    ) -> List[Material]:
        if self._GeneratedItemType == Material:
            return items
        return [Material.create(self._convert_generated_item(item)) for item in items]

    @staticmethod
    def _convert_generated_item(item: _GeneratedItemType):
        material_config = item
        return material_config

    def _finalize(self, materials: List[Material], configuration: _ConfigurationType) -> List[Material]:
        materials_with_metadata = [self._update_material_metadata(material, configuration) for material in materials]
        return [self._update_material_name(material, configuration) for material in materials_with_metadata]

    def get_materials(
        self,
        configuration: _ConfigurationType,
        selector_parameters: Optional[_SelectorParametersType] = None,
        post_process_parameters: Optional[_PostProcessParametersType] = None,
    ) -> List[Material]:
        generated_items = self._generate_or_get_from_cache(configuration)
        sorted_items = self._sort(generated_items)
        selected_items = self._select(sorted_items, selector_parameters)
        materials = self._post_process(selected_items, post_process_parameters)
        finalized_materials = self._finalize(materials, configuration)
        return finalized_materials

    def get_material(
        self,
        configuration: _ConfigurationType,
        selector_parameters: Optional[_SelectorParametersType] = None,
        post_process_parameters: Optional[_PostProcessParametersType] = None,
    ) -> Material:
        return self.get_materials(configuration, selector_parameters, post_process_parameters)[
            self.selector_parameters.default_index
        ]

    def _update_material_name(self, material, configuration) -> Material:
        # Do nothing by default
        return material

    def _update_material_metadata(self, material, configuration) -> Material:
        if "build" not in material.metadata:
            material.metadata["build"] = {}
        material.metadata["build"]["configuration"] = configuration.to_json()
        return material
