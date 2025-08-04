from pydantic import Field

from mat3ra.made.material import Material

from .material_build_metadata import MaterialBuildMetadata


class MaterialWithBuildMetadata(Material):
    metadata: MaterialBuildMetadata = Field(default_factory=MaterialBuildMetadata)
