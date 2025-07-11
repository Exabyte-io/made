from typing import List, Optional

from pydantic import Field

from mat3ra.made.metadata import BaseMetadata


class BuildMetadata(BaseMetadata):
    configuration: dict = Field(default_factory=dict)
    build_parameters: dict = Field(default_factory=dict)


class MaterialMetadata(BaseMetadata):
    # Other metadata fields can be added as needed
    build: Optional[List[BuildMetadata]] = Field([])
