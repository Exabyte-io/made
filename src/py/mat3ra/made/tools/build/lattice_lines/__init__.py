from typing import Optional, Tuple, Union

from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.lattice_lines import CrystalLatticeLinesMaterialAnalyzer
from .builders import (
    CrystalLatticeLinesBuilder,
    CrystalLatticeLinesRepeatedBuilder,
)
from .configuration import (
    CrystalLatticeLinesConfiguration,
    CrystalLatticeLinesUniqueRepeatedConfiguration,
    get_miller_indices_from_edge_type,
    EdgeTypes,
)
from .. import MaterialWithBuildMetadata
from ..slab.termination_utils import select_slab_termination


def create_lattice_lines_config_and_material(
    material: Union[Material, MaterialWithBuildMetadata],
    miller_indices_2d: Optional[Tuple[int, int]],
    edge_type: Optional[EdgeTypes],
    width: int,
    length: int,
    termination_formula: Optional[str] = None,
):
    if miller_indices_2d is None and edge_type is None:
        raise ValueError("Either miller_indices_2d or edge_type must be provided")
    if miller_indices_2d is None and edge_type is not None:
        miller_indices_2d = get_miller_indices_from_edge_type(edge_type)

    lattice_lines_analyzer = CrystalLatticeLinesMaterialAnalyzer(material=material, miller_indices_2d=miller_indices_2d)
    terminations = lattice_lines_analyzer.terminations
    termination = select_slab_termination(terminations, termination_formula)

    lattice_lines_config = CrystalLatticeLinesUniqueRepeatedConfiguration(
        crystal=material,
        miller_indices_2d=miller_indices_2d,
        termination_top=termination,
        number_of_repetitions_width=width,
        number_of_repetitions_length=length,
    )
    return lattice_lines_config


__all__ = [
    "CrystalLatticeLinesConfiguration",
    "CrystalLatticeLinesUniqueRepeatedConfiguration",
    "CrystalLatticeLinesBuilder",
    "CrystalLatticeLinesRepeatedBuilder",
    "create_lattice_lines_config_and_material",
]
