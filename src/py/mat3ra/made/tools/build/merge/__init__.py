# Compatibility layer for build.merge
from ...build_components.operations.core.combinations.merge.builder import MergeBuilder
from ...build_components.operations.core.combinations.merge.configuration import MergeConfiguration
from ...build_components.entities.reusable.three_dimensional.crystal_lattice_base.build_parameters import BaseBuilderParameters as MergeBuilderParameters

__all__ = [
    "MergeBuilder",
    "MergeConfiguration",
    "MergeBuilderParameters",
]