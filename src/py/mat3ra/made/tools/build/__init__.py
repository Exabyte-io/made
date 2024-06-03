from typing import List, Optional, Any

from pydantic import BaseModel

from ...material import Material


class BaseBuilder(BaseModel):
    build_parameters: Any = None
    _BuildParametersType: Any = None
    _DefaultBuildParameters: Any = None

    _ConfigurationType: Any = Any
    _GeneratedItemType: Any = Any
    _SelectorParametersType: Any = None
    _PostProcessParametersType: Any = None

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
        return [Material(self._convert_generated_item(item)) for item in items]

    @staticmethod
    def _convert_generated_item(item: _GeneratedItemType):
        material_config = item
        return material_config

    def _finalize(self, materials: List[Material], configuration: _ConfigurationType) -> List[Material]:
        return [self._update_material_name(material, configuration) for material in materials]

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
        return self.get_materials(configuration, selector_parameters, post_process_parameters)[0]

    def _update_material_name(self, material, configuration):
        # Do nothing by default
        return material
