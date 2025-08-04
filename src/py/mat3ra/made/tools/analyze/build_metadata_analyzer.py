from __future__ import annotations

from typing import Generic, Type, TypeVar

from ..build import BuildMetadata, MaterialWithBuildMetadata
from . import BaseMaterialAnalyzer

TConfiguration = TypeVar("TConfiguration")
TBuildParameters = TypeVar("TBuildParameters")


class BuildMetadataAnalyzer(BaseMaterialAnalyzer, Generic[TConfiguration, TBuildParameters]):
    material: MaterialWithBuildMetadata
    configuration_cls: Type[TConfiguration]
    build_parameters_cls: Type[TBuildParameters]

    @property
    def build_metadata(self) -> BuildMetadata:
        metadata = self.material.metadata.get_build_metadata_of_type(self.configuration_cls.__name__)
        if not metadata:
            raise ValueError(f"Build metadata for {self.configuration_cls.__name__} not found.")
        return metadata

    @property
    def build_configuration(self) -> TConfiguration:
        return self.configuration_cls(**self.build_metadata.configuration)

    @property
    def build_parameters(self) -> TBuildParameters:
        return self.build_parameters_cls(**self.build_metadata.build_parameters)
