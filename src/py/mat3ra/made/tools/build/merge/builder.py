# Compatibility layer
from ...build_components.operations.core.combinations.merge.builder import MergeBuilder
from mat3ra.made.tools.build_components.entities.reusable.base_builder import (
    BaseBuilderParameters as MergeBuilderParameters,
)

__all__ = ["MergeBuilder", "MergeBuilderParameters"]
