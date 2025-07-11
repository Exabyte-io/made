from typing import List, Optional

from pydantic import Field

from mat3ra.made.material import Material
from mat3ra.made.metadata import BaseMetadata


class BuildMetadata(BaseMetadata):
    configuration: dict = Field(default_factory=dict)
    build_parameters: dict = Field(default_factory=dict)

    @staticmethod
    def set_build_metadata(
        material: Material,
        configuration: dict,
        build_parameters: Optional[dict] = None,
    ) -> None:
        metadata = MaterialMetadata(**material.metadata or {})
        build_metadata = BuildMetadata(configuration=configuration, build_parameters=build_parameters or {})
        metadata.build.append(build_metadata)
        material.metadata = metadata.to_dict()


class MaterialMetadata(BaseMetadata):
    # Other metadata fields can be added as needed
    build: Optional[List[BuildMetadata]] = Field([])
