import json
from typing import Any, Dict

from pydantic import Field

from mat3ra.code.entity import InMemoryEntityPydantic


class BaseMetadata(InMemoryEntityPydantic):
    model_config = {"arbitrary_types_allowed": True, "extra": "allow"}

    @staticmethod
    def _to_dict(obj: Any) -> dict:
        if not obj:
            return {}

        if hasattr(obj, "to_dict") and callable(obj.to_dict):
            return obj.to_dict()
        if hasattr(obj, "to_json") and callable(obj.to_json):
            data = obj.to_json()
            if isinstance(data, str):
                return json.loads(data)
            return data
        return {}

    def update(self, **kwargs: Any):
        for key, value in kwargs.items():
            if hasattr(self, key) and isinstance(getattr(self, key), dict):
                getattr(self, key).update(self._to_dict(value))


class BuildMetadata(BaseMetadata):
    configuration: Dict[str, Any] = Field(default_factory=dict)
    build_parameters: Dict[str, Any] = Field(default_factory=dict)


class MaterialMetadata(BaseMetadata):
    # Other metadata fields can be added as needed
    build: BuildMetadata = Field(default_factory=BuildMetadata)
