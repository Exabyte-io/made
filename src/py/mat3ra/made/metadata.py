import json
from typing import Any

from mat3ra.code.entity import InMemoryEntityPydantic


def to_dict(obj: Any) -> dict:
    if isinstance(obj, dict):
        return obj
    if hasattr(obj, "to_dict") and callable(obj.to_dict):
        return obj.to_dict()
    # TODO: remove conditional checks below when all configurations moved to Pydantic
    if hasattr(obj, "to_json") and callable(obj.to_json):
        data = obj.to_json()
        if isinstance(data, str):
            return json.loads(data)
        return data
    return {}


class BaseMetadata(InMemoryEntityPydantic):
    model_config = {"arbitrary_types_allowed": True, "extra": "allow"}

    def update(self, **kwargs: Any):
        for key, value in kwargs.items():
            if hasattr(self, key) and isinstance(getattr(self, key), dict):
                getattr(self, key).update(to_dict(value))
