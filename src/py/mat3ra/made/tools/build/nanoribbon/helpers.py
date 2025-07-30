from typing import Tuple, Optional, TypeVar, Union

from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum

from mat3ra.made.material import Material
from mat3ra.made.tools.build import BaseBuilder, BaseBuilderParameters, MaterialWithBuildMetadata
from . import NanoribbonConfiguration
from .builders import NanoribbonBuilder, NanoribbonBuilderParameters
from ..lattice_lines import create_lattice_lines_config_and_material
from ..lattice_lines.configuration import EdgeTypes
from ..nanotape.configuration import NanoTapeConfiguration
from ..vacuum.configuration import VacuumConfiguration

T = TypeVar("T", bound=BaseBuilder)
P = TypeVar("P", bound=BaseBuilderParameters)


def create_nanoribbon(
    material: Union[Material, MaterialWithBuildMetadata],
    miller_indices_2d: Optional[Tuple[int, int]] = None,
    edge_type: EdgeTypes = EdgeTypes.zigzag,
    width: int = 2,
    length: int = 2,
    vacuum_width: float = 10.0,
    vacuum_length: float = 10.0,
    use_rectangular_lattice: bool = True,
    termination_formula: Optional[str] = None,
) -> Material:
    """
    Creates a nanoribbon material from a monolayer material.

    Args:
        material: The monolayer material to create the nanoribbon from (assumes vacuum is present).
        miller_indices_2d: The (u,v) Miller indices for the nanoribbon direction.
        edge_type: Edge type string ("zigzag"/"armchair"). Optional if miller_indices_2d is provided.
        width: The width of the nanoribbon in number of unit cells.
        length: The length of the nanoribbon in number of unit cells.
        vacuum_width: The width of the vacuum region in Angstroms (cartesian).
        vacuum_length: The length of the vacuum region in Angstroms (cartesian).
        use_rectangular_lattice: Whether the nanoribbon is rectangular.
        termination_formula: The termination formula to use for the nanoribbon (e.g., "Si").

    Returns:
        Material: The generated nanoribbon material.
    """
    lattice_lines_config = create_lattice_lines_config_and_material(
        material=material,
        miller_indices_2d=miller_indices_2d,
        edge_type=edge_type,
        width=width,
        length=length,
        termination_formula=termination_formula,
    )
    nanotape_vacuum_config = VacuumConfiguration(
        size=vacuum_width,
        direction=AxisEnum.y,
    )
    nanotape_config = NanoTapeConfiguration(
        stack_components=[lattice_lines_config, nanotape_vacuum_config],
        direction=AxisEnum.y,
    )
    vacuum_config = VacuumConfiguration(
        size=vacuum_length,
        direction=AxisEnum.x,
    )
    config = NanoribbonConfiguration(
        stack_components=[nanotape_config, vacuum_config],
        direction=AxisEnum.x,
    )
    builder = NanoribbonBuilder(
        build_parameters=NanoribbonBuilderParameters(use_rectangular_lattice=use_rectangular_lattice)
    )
    return builder.get_material(config)
