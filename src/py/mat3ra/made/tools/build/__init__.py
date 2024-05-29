from typing import List, Any

from ...material import Material


class BaseBuilder:
    def __generator(self, **kwargs) -> Any:
        raise NotImplementedError

    def get_materials(
        self,
        configuration,
        **kwargs,
    ) -> List[Material]:
        raise NotImplementedError

    def get_material(self, configuration, **kwargs) -> Material:
        raise NotImplementedError
