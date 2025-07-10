from typing import List

from pydantic import Field, model_validator

from mat3ra.made.metadata import BaseMetadata


class BuildMetadata(BaseMetadata):
    configuration: dict = Field(default_factory=dict)
    build_parameters: dict = Field(default_factory=dict)


class MaterialMetadata(BaseMetadata):
    # Other metadata fields can be added as needed
    build: List[BuildMetadata] = Field(default_factory=lambda: [BuildMetadata()])
