from typing import Any, Dict, Type, TypeVar

from pydantic import BaseModel, Field

MetadataChild = TypeVar("MetadataChild", bound="BaseMetadata")


class BaseMetadata(BaseModel):
    model_config = {"extra": "allow"}

    @classmethod
    def from_dict(cls: Type[MetadataChild], data: dict) -> MetadataChild:
        return cls.model_validate(data or {})


class BuildMetadata(BaseMetadata):
    configuration: Dict[str, Any] = Field(default_factory=dict)
    build_parameters: Dict[str, Any] = Field(default_factory=dict)

    def update(self, configuration: Any, build_parameters: Any = None):
        # TODO: remove conditional checks when all configurations moved to Pydantic
        if hasattr(configuration, "to_dict") and callable(configuration.to_dict):
            config_dict = configuration.to_dict()
        elif hasattr(configuration, "to_json") and callable(configuration.to_json):
            config_dict = configuration.to_json()
        else:
            config_dict = {}

        self.configuration.update(config_dict)

        build_parameters_dict = {}
        if hasattr(build_parameters, "to_dict"):
            build_parameters_dict = build_parameters.model_dump()
        elif hasattr(build_parameters, "to_json"):
            build_parameters_dict = build_parameters.to_json()

        self.build_parameters.update(build_parameters_dict)


class MaterialMetadata(BaseMetadata):
    # Other metadata fields can be added as needed
    build: BuildMetadata = Field(default_factory=BuildMetadata)
