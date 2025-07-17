import json
from typing import Any

from mat3ra.code.entity import InMemoryEntityPydantic
from pydantic import model_validator


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

    @model_validator(mode="before")
    def convert_fields_to_dict(cls, values):
        for key, value in values.items():
            field = cls.model_fields.get(key)
            if field and field.annotation is dict and value is None:
                values[key] = {}
            elif field and field.annotation is dict and value is not None and not isinstance(value, dict):
                values[key] = to_dict(value)
            elif isinstance(value, InMemoryEntityPydantic):
                values[key] = to_dict(value)
            elif hasattr(value, "to_dict") and callable(value.to_dict):
                values[key] = to_dict(value)
            # TODO: remove when Pydantic configurations are fully implemented
            elif hasattr(value, "to_json") and callable(value.to_json):
                values[key] = to_dict(value)
        return values

    def update(self, **kwargs: Any):
        for key, value in kwargs.items():
            if hasattr(self, key) and isinstance(getattr(self, key), dict):
                getattr(self, key).update(to_dict(value))
