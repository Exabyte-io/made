from typing import Optional, Tuple
from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.lattice_lines import CrystalLatticeLinesMaterialAnalyzer
from .lattice_lines.configuration import (
    CrystalLatticeLinesUniqueRepeatedConfiguration,
    get_miller_indices_from_edge_type,
    EdgeTypes,
)
from .lattice_lines.builders import CrystalLatticeLinesRepeatedBuilder
from .vacuum.configuration import VacuumConfiguration
from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from mat3ra.made.tools.build.slab.entities import Termination

def create_lattice_lines_config_and_material(
    material: Material,
    miller_indices_2d: Optional[Tuple[int, int]],
    edge_type: Optional[EdgeTypes],
    width: int,
    length: int,
    termination: Optional[Termination] = None,
):
    if miller_indices_2d is None and edge_type is None:
        raise ValueError("Either miller_indices_2d or edge_type must be provided")
    if miller_indices_2d is None and edge_type is not None:
        miller_indices_2d = get_miller_indices_from_edge_type(edge_type)

    lattice_lines_analyzer = CrystalLatticeLinesMaterialAnalyzer(
        material=material, miller_indices_2d=miller_indices_2d
    )
    if termination is None:
        termination = lattice_lines_analyzer.default_termination

    lattice_lines_config = CrystalLatticeLinesUniqueRepeatedConfiguration(
        crystal=material,
        miller_indices_2d=miller_indices_2d,
        termination_top=termination,
        number_of_repetitions_width=width,
        number_of_repetitions_length=length,
    )
    lattice_lines_builder = CrystalLatticeLinesRepeatedBuilder()
    lattice_lines_material = lattice_lines_builder.get_material(lattice_lines_config)
    return lattice_lines_config, lattice_lines_material

def create_vacuum_config(
    size: float,
    crystal: Material,
    direction: AxisEnum,
) -> VacuumConfiguration:
    return VacuumConfiguration(size=size, crystal=crystal, direction=direction) 