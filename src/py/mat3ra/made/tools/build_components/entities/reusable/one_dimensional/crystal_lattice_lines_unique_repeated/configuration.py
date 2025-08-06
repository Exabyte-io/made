from typing import Optional

from ....auxiliary.two_dimensional.termination import Termination
from ..crystal_lattice_lines.configuration import CrystalLatticeLinesConfiguration


class CrystalLatticeLinesUniqueRepeatedConfiguration(CrystalLatticeLinesConfiguration):
    """
    Configuration for creating repeated crystal lattice lines with termination.

    Args:
        crystal: The monolayer material to create the lattice lines from.
        miller_indices_2d: The (u,v) Miller indices for the line direction.
        termination_top: The termination to use for the lattice lines.
    """

    # TODO: right and left for x terminations
    termination_top: Termination
    termination_bottom: Optional[Termination] = None
    number_of_repetitions_width: int = 1
    number_of_repetitions_length: int = 1
