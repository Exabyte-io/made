# Compatibility layer for build.merge
from ...build_components.operations.core.combinations.merge.builder import MergeBuilder
from ...build_components.operations.core.combinations.merge.configuration import MergeConfiguration
from mat3ra.made.tools.build_components.entities.reusable.base_builder import (
    BaseBuilderParameters as MergeBuilderParameters,
)

__all__ = [
    "MergeBuilder",
    "MergeConfiguration",
    "MergeBuilderParameters",
]
