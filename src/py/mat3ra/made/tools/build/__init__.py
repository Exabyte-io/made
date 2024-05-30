from typing import List, Any, Optional, Dict

from ...material import Material


class BaseBuilder:
    __ConfigurationType: Any = Any
    __GeneratedItemType: Any = Any
    __BuildParametersType: Optional[Any] = Optional[Any]
    __SelectorParametersType: Optional[Dict] = Optional[Dict]
    __PostProcessParametersType: Optional[Dict] = Optional[Dict]

    def __init__(self, build_parameters: __BuildParametersType = None) -> None:
        self.build_parameters = build_parameters
        self.__generated_items: List[List[Any]] = []
        self.__configurations: List[Any] = []

    def __generate_or_get_from_cache(self, configuration: __ConfigurationType) -> List[Any]:
        if configuration not in self.__configurations:
            self.__configurations.append(configuration)
            self.__generated_items.append(self.__generate(configuration))
        return self.__generated_items[self.__configurations.index(configuration)]

    def __generate(self, configuration: __ConfigurationType) -> List[__GeneratedItemType]:
        return []

    def __sort(self, items: List[__GeneratedItemType]) -> List[__GeneratedItemType]:
        return items

    def __select(
        self, items: List[__GeneratedItemType], selector_parameters: __SelectorParametersType
    ) -> List[__GeneratedItemType]:
        return items

    def __post_process(
        self, items: List[__GeneratedItemType], post_process_parameters: __PostProcessParametersType
    ) -> List[Material]:
        return items

    def get_materials(
        self,
        configuration: __ConfigurationType,
        selector_parameters: __SelectorParametersType = None,
        post_process_parameters: __PostProcessParametersType = None,
    ) -> List[Material]:
        generated_items = self.__generate_or_get_from_cache(configuration)
        sorted_items = self.__sort(generated_items)
        selected_items = self.__select(sorted_items, selector_parameters)
        return self.__post_process(selected_items, post_process_parameters)

    def get_material(
        self,
        configuration: __ConfigurationType,
        selector_parameters: __SelectorParametersType = None,
        post_process_parameters: __PostProcessParametersType = None,
    ) -> Material:
        return self.get_materials(configuration, selector_parameters, post_process_parameters)[0]
