from typing import List, Optional

from pydantic import Field

from mat3ra.made.material import Material
from mat3ra.made.metadata import BaseMetadata


class BuildMetadata(BaseMetadata):
    configuration: dict = Field(default_factory=dict)
    build_parameters: dict = Field(default_factory=dict)


class MaterialBuildMetadata(BaseMetadata):
    # Other metadata fields can be added as needed
    build: List[BuildMetadata] = Field([])

    def get_build_metadata_of_type(self, configuration_type: str) -> Optional[BuildMetadata]:
        """
        Returns the first build configuration of the specified type.
        Args:
            configuration_type (str): The type of configuration to search for.
        Returns:
            Optional[dict]: The build configuration if found, otherwise None.
        """
        for build_step in reversed(self.build):
            if build_step.configuration.get("type") == configuration_type:
                return BuildMetadata(
                    configuration=build_step.configuration, build_parameters=build_step.build_parameters
                )
        return None

    def add_build_metadata_step(self, configuration: dict, build_parameters: Optional[dict] = None) -> None:
        build_metadata = BuildMetadata(configuration=configuration, build_parameters=build_parameters)
        self.build.append(build_metadata)


class MaterialWithBuildMetadata(Material):
    metadata: MaterialBuildMetadata = Field(default_factory=MaterialBuildMetadata)
