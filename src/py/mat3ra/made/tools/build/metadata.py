from typing import List

from pydantic import Field

from mat3ra.made.metadata import BaseMetadata


class BuildMetadata(BaseMetadata):
    configuration: dict = Field(default_factory=dict)
    build_parameters: dict = Field(default_factory=dict)


class MaterialMetadata(BaseMetadata):
    # Other metadata fields can be added as needed
    build: List[BuildMetadata] = Field(default_factory=lambda: [BuildMetadata()])


def get_slab_build_configuration(metadata: dict):
    """
    Extract slab build configuration from material metadata.

    Args:
        metadata: Material metadata dictionary

    Returns:
        SlabConfiguration: The slab configuration from the material's metadata

    Raises:
        ValueError: If the material is not a slab
    """
    from .slab.configurations import SlabConfiguration  # Import here to avoid circular imports

    material_metadata = MaterialMetadata(**metadata)
    slab_build_configuration_dict = material_metadata.build[-1].configuration

    if slab_build_configuration_dict.get("type") != "SlabConfiguration":
        raise ValueError("Material is not a slab.")

    return SlabConfiguration(**slab_build_configuration_dict)
