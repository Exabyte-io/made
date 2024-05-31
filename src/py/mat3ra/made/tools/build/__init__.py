from typing import List, Optional, Any

from ...material import Material


class BaseBuilder:
    _ConfigurationType: Any = Any
    _BuildParametersType: Any = Any
    _GeneratedItemType: Any = Any
    _SelectorParametersType: Any = Any
    _PostProcessParametersType: Any = Any

    def __init__(self, build_parameters: _BuildParametersType) -> None:
        self.build_parameters = build_parameters
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
        self, items: List[_GeneratedItemType], selector_parameters: Optional[_SelectorParametersType] = None
    ) -> List[_GeneratedItemType]:
        return items

    def _post_process(
        self, items: List[_GeneratedItemType], post_process_parameters: Optional[_PostProcessParametersType] = None
    ) -> List[Material]:
        # TODO: convert from _GeneratedItemType to Material
        return [Material(item) for item in items]

    def _update_material_name(self, material: Material, configuration: _ConfigurationType) -> Material:
        return material

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
        materials_with_name = [self._update_material_name(material, configuration) for material in materials]
        return materials_with_name

    def get_material(
        self,
        configuration: _ConfigurationType,
        selector_parameters: Optional[_SelectorParametersType] = None,
        post_process_parameters: Optional[_PostProcessParametersType] = None,
    ) -> Material:
        return self.get_materials(configuration, selector_parameters, post_process_parameters)[0]
