from mat3ra.made.metadata import BaseMetadata
from pydantic import Field


class BuildMetadata(BaseMetadata):
    configuration: dict = Field(default_factory=dict)
    build_parameters: dict = Field(default_factory=dict)
