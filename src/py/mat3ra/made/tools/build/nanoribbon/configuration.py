from typing import Optional, Tuple

from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum

from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.lattice_lines import CrystalLatticeLinesAnalyzer
from mat3ra.made.tools.build.nanoribbon.builders import (
    CrystalLatticeLinesRepeatedBuilder,
    NanoTapeBuilder,
)
from .crystal_lattice_lines_configuration import (
    CrystalLatticeLinesUniqueRepeatedConfiguration,
)
from .enums import EdgeTypes
from .nano_tape_configuration import NanoTapeConfiguration
from .nanoribbon_configuration import NanoribbonConfiguration
from ..slab.entities import Termination
from ..vacuum.configuration import VacuumConfiguration


# Helper functions for shorthand notation
def get_miller_indices_from_edge_type(edge_type: EdgeTypes) -> Tuple[int, int]:
    """
    Convert edge type shorthand to (u,v) Miller indices.

    Args:
        edge_type: "zigzag" or "armchair"

    Returns:
        Tuple of (u,v) Miller indices.
    """
    if edge_type.lower() == EdgeTypes.zigzag:
        return (1, 1)
    elif edge_type.lower() == EdgeTypes.armchair:
        return (0, 1)
    else:
        raise ValueError(f"Unknown edge type: {edge_type}. Use 'zigzag' or 'armchair'.")


# Factory methods for configurations
def create_nanotape_configuration(
    material: Material,
    miller_indices_uv: Optional[Tuple[int, int]] = None,
    edge_type: Optional[EdgeTypes] = EdgeTypes.zigzag,
    width: int = 2,
    length: int = 2,
    vacuum_width: float = 10.0,
    termination: Optional[Termination] = None,
) -> NanoTapeConfiguration:
    if miller_indices_uv is None and edge_type is None:
        raise ValueError("Either miller_indices_uv or edge_type must be provided")
    if miller_indices_uv is None and edge_type is not None:
        miller_indices_uv = get_miller_indices_from_edge_type(edge_type)

    lattice_lines_analyzer = CrystalLatticeLinesAnalyzer(material=material, miller_indices_uv=miller_indices_uv)
    if termination is None:
        termination = lattice_lines_analyzer.default_termination

    lattice_lines_config = CrystalLatticeLinesUniqueRepeatedConfiguration(
        crystal=material,
        miller_indices_uv=miller_indices_uv,
        termination_top=termination,
        number_of_repetitions_width=width,
        number_of_repetitions_length=length,
    )
    lattice_lines_builder = CrystalLatticeLinesRepeatedBuilder()
    lattice_lines_material = lattice_lines_builder.get_material(lattice_lines_config)
    vacuum_config = VacuumConfiguration(
        size=vacuum_width,
        crystal=lattice_lines_material,
        direction=AxisEnum.y,
    )
    return NanoTapeConfiguration(
        stack_components=[lattice_lines_config, vacuum_config],
        direction=AxisEnum.y,
    )


def create_nanoribbon_configuration(
    material: Material,
    miller_indices_uv: Optional[Tuple[int, int]] = None,
    edge_type: Optional[EdgeTypes] = EdgeTypes.zigzag,
    width: int = 2,
    length: int = 2,
    vacuum_width: float = 10.0,
    vacuum_length: float = 10.0,
    termination: Optional[Termination] = None,
) -> NanoribbonConfiguration:
    if miller_indices_uv is None and edge_type is None:
        raise ValueError("Either miller_indices_uv or edge_type must be provided")
    if miller_indices_uv is None and edge_type is not None:
        miller_indices_uv = get_miller_indices_from_edge_type(edge_type)

    nanotape_config = create_nanotape_configuration(
        material=material,
        miller_indices_uv=miller_indices_uv,
        width=width,
        length=length,
        termination=termination,
        vacuum_width=vacuum_width,
    )
    nanotape_builder = NanoTapeBuilder()
    nanotape_material = nanotape_builder.get_material(nanotape_config)
    vacuum_config = VacuumConfiguration(
        size=vacuum_length,
        crystal=nanotape_material,
        direction=AxisEnum.x,
    )
    return NanoribbonConfiguration(
        stack_components=[nanotape_config, vacuum_config],
        direction=AxisEnum.x,
    )
