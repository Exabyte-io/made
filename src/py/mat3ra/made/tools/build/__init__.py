from typing import List, Optional, TypeVar, Generic

from ...material import Material


ConfigurationType = TypeVar("ConfigurationType")
BuildParametersType = TypeVar("BuildParametersType")
GeneratedItemType = TypeVar("GeneratedItemType")
SelectorParametersType = TypeVar("SelectorParametersType")
PostProcessParametersType = TypeVar("PostProcessParametersType")


class BaseBuilder(
    Generic[
        ConfigurationType, BuildParametersType, GeneratedItemType, SelectorParametersType, PostProcessParametersType
    ]
):

    def __init__(self, build_parameters: BuildParametersType) -> None:
        self.build_parameters = build_parameters
        self.__generated_items: List[List[GeneratedItemType]] = []
        self.__configurations: List[ConfigurationType] = []

    def __generate(self, configuration: ConfigurationType) -> List[GeneratedItemType]:
        return []

    def __generate_or_get_from_cache(self, configuration: ConfigurationType) -> List[GeneratedItemType]:
        if configuration not in self.__configurations:
            self.__configurations.append(configuration)
            self.__generated_items.append(self.__generate(configuration))
        return self.__generated_items[self.__configurations.index(configuration)]

    def __sort(self, items: List[GeneratedItemType]) -> List[GeneratedItemType]:
        return items

    def __select(
        self, items: List[GeneratedItemType], selector_parameters: Optional[SelectorParametersType] = None
    ) -> List[GeneratedItemType]:
        return items

    def __post_process(
        self, items: List[GeneratedItemType], post_process_parameters: Optional[PostProcessParametersType] = None
    ) -> List[Material]:
        return [Material(item) for item in items]

    def get_materials(
        self,
        configuration: ConfigurationType,
        selector_parameters: Optional[SelectorParametersType] = None,
        post_process_parameters: Optional[PostProcessParametersType] = None,
    ) -> List[Material]:
        generated_items = self.__generate_or_get_from_cache(configuration)
        sorted_items = self.__sort(generated_items)
        selected_items = self.__select(sorted_items, selector_parameters)
        return self.__post_process(selected_items, post_process_parameters)

    def get_material(
        self,
        configuration: ConfigurationType,
        selector_parameters: Optional[SelectorParametersType] = None,
        post_process_parameters: Optional[PostProcessParametersType] = None,
    ) -> Material:
        return self.get_materials(configuration, selector_parameters, post_process_parameters)[0]
