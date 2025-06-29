from typing import Tuple, Optional
from mat3ra.made.material import Material
from mat3ra.made.tools.build import BaseConfigurationPydantic
from ..slab.entities import Termination


class CrystalLatticeLinesConfiguration(BaseConfigurationPydantic):
    """
    Configuration for creating crystal lattice lines from a material.

    Args:
        crystal: The monolayer material to create the lattice lines from.
        miller_indices_uv: The (u,v) Miller indices for the line direction.
    """

    crystal: Material
    miller_indices_uv: Tuple[int, int]


class CrystalLatticeLinesUniqueRepeatedConfiguration(CrystalLatticeLinesConfiguration):
    """
    Configuration for creating repeated crystal lattice lines with termination.

    Args:
        crystal: The monolayer material to create the lattice lines from.
        miller_indices_uv: The (u,v) Miller indices for the line direction.
        termination_top: The termination to use for the lattice lines.
    """

    # TODO: right and left for x terminations
    termination_top: Termination
    termination_bottom: Optional[Termination] = None
    number_of_repetitions_width: int = 1
    number_of_repetitions_length: int = 1 