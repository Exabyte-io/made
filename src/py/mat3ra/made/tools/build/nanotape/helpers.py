from typing import Tuple, Optional, Union

from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum

from mat3ra.made.material import Material
from . import NanoTapeConfiguration
from .builders import NanoTapeBuilder, NanoTapeBuilderParameters
from .. import MaterialWithBuildMetadata
from ..lattice_lines.configuration import EdgeTypes
from ..lattice_lines import create_lattice_lines_config_and_material
from ..vacuum.configuration import VacuumConfiguration


def create_nanotape(
    material: Union[Material, MaterialWithBuildMetadata],
    miller_indices_2d: Optional[Tuple[int, int]] = None,
    edge_type: EdgeTypes = EdgeTypes.zigzag,
    width: int = 2,
    length: int = 2,
    vacuum_width: float = 10.0,
    use_rectangular_lattice: bool = True,
    termination_formula: Optional[str] = None,
) -> Material:
    """
    Creates a nanotape material from a monolayer material.

    Args:
        material: The monolayer material to create the nanotape from (assumes vacuum is present).
        miller_indices_2d: The (u,v) Miller indices for the nanotape direction.
        edge_type: Edge type string ("zigzag"/"armchair"). Optional if miller_indices_2d is provided.
        width: The width of the nanotape in number of unit cells.
        length: The length of the nanotape in number of unit cells.
        vacuum_width: The width of the vacuum region in Angstroms (cartesian).
        use_rectangular_lattice: Whether the nanotape is rectangular.
        termination_formula: The termination formula to use for the nanotape (e.g., "Si").

    Returns:
        Material: The generated nanotape material.
    """
    lattice_lines_config, lattice_lines_material = create_lattice_lines_config_and_material(
        material=material,
        miller_indices_2d=miller_indices_2d,
        edge_type=edge_type,
        width=width,
        length=length,
        termination_formula=termination_formula,
    )
    vacuum_config = VacuumConfiguration(
        size=vacuum_width,
        crystal=lattice_lines_material,
        direction=AxisEnum.y,
    )
    config = NanoTapeConfiguration(
        stack_components=[lattice_lines_config, vacuum_config],
        direction=AxisEnum.y,
    )
    builder = NanoTapeBuilder(
        build_parameters=NanoTapeBuilderParameters(use_rectangular_lattice=use_rectangular_lattice)
    )
    return builder.get_material(config)
