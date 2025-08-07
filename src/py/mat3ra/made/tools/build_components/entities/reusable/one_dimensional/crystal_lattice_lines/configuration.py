from typing import Tuple, Union

from mat3ra.made.material import Material
from mat3ra.made.tools.build_components.entities.reusable.base_builder import BaseConfigurationPydantic

from ..... import MaterialWithBuildMetadata


class CrystalLatticeLinesConfiguration(BaseConfigurationPydantic):
    """
    Configuration for creating crystal lattice lines from a material.

    Args:
        crystal: The monolayer material to create the lattice lines from.
        miller_indices_2d: The (u,v) Miller indices for the line direction.
    """

    crystal: Union[Material, MaterialWithBuildMetadata]
    miller_indices_2d: Tuple[int, int]
