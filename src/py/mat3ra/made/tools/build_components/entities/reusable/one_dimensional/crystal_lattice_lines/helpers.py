from typing import Optional, Tuple, Union

from mat3ra.made.material import Material

from ......analyze.lattice_lines import CrystalLatticeLinesMaterialAnalyzer
from ......build.pristine_structures.two_dimensional.slab.termination_utils import select_slab_termination
from ..... import MaterialWithBuildMetadata
from ..crystal_lattice_lines_unique_repeated.configuration import CrystalLatticeLinesUniqueRepeatedConfiguration
from .edge_types import EdgeTypesEnum, get_miller_indices_from_edge_type


def create_lattice_lines_config_and_material(
    material: Union[Material, MaterialWithBuildMetadata],
    miller_indices_2d: Optional[Tuple[int, int]],
    edge_type: Optional[EdgeTypesEnum],
    width: int,
    length: int,
    termination_formula: Optional[str] = None,
):
    """
    Create a lattice lines configuration from a material.

    Args:
        material: The monolayer material to create the lattice lines from.
        miller_indices_2d: The (u,v) Miller indices for the line direction.
        edge_type: Edge type ("zigzag" or "armchair") - alternative to miller_indices_2d.
        width: Number of repetitions in width direction.
        length: Number of repetitions in length direction.
        termination_formula: Formula of the termination to use.

    Returns:
        CrystalLatticeLinesUniqueRepeatedConfiguration: The configuration object.
    """
    if miller_indices_2d is None and edge_type is None:
        raise ValueError("Either miller_indices_2d or edge_type must be provided")

    resolved_miller_indices: Tuple[int, int]
    if miller_indices_2d is not None:
        resolved_miller_indices = miller_indices_2d
    elif edge_type is not None:
        resolved_miller_indices = get_miller_indices_from_edge_type(edge_type)
    else:
        raise ValueError("Either miller_indices_2d or edge_type must be provided")

    lattice_lines_analyzer = CrystalLatticeLinesMaterialAnalyzer(
        material=material, miller_indices_2d=resolved_miller_indices
    )
    terminations = lattice_lines_analyzer.terminations
    termination = select_slab_termination(terminations, termination_formula)

    lattice_lines_config = CrystalLatticeLinesUniqueRepeatedConfiguration(
        crystal=material,
        miller_indices_2d=resolved_miller_indices,
        termination_top=termination,
        number_of_repetitions_width=width,
        number_of_repetitions_length=length,
    )
    return lattice_lines_config
