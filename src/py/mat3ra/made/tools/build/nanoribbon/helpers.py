from typing import Tuple, Optional, TypeVar

from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum

from mat3ra.made.material import Material
from mat3ra.made.tools.build import BaseBuilder, BaseBuilderParameters
from mat3ra.made.tools.build.slab.entities import Termination
from . import NanoribbonConfiguration
from .builders import NanoribbonBuilder, NanoribbonBuilderParameters
from ..lattice_lines.configuration import EdgeTypes
from ..lattice_lines_utils import create_lattice_lines_config_and_material, create_vacuum_config
from ..nanotape.builders import NanoTapeBuilder
from ..nanotape.configuration import NanoTapeConfiguration

T = TypeVar("T", bound=BaseBuilder)
P = TypeVar("P", bound=BaseBuilderParameters)


def create_nanoribbon(
    material: Material,
    miller_indices_2d: Optional[Tuple[int, int]] = None,
    edge_type: EdgeTypes = EdgeTypes.zigzag,
    width: int = 2,
    length: int = 2,
    vacuum_width: float = 10.0,
    vacuum_length: float = 10.0,
    use_rectangular_lattice: bool = True,
    termination: Optional[Termination] = None,
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
        termination: The termination to use for the nanoribbon. If None, uses default termination.

    Returns:
        Material: The generated nanoribbon material.
    """
    lattice_lines_config, lattice_lines_material = create_lattice_lines_config_and_material(
        material=material,
        miller_indices_2d=miller_indices_2d,
        edge_type=edge_type,
        width=width,
        length=length,
        termination=termination,
    )
    nanotape_vacuum_config = create_vacuum_config(
        size=vacuum_width,
        crystal=lattice_lines_material,
        direction=AxisEnum.y,
    )
    nanotape_config = NanoTapeConfiguration(
        stack_components=[lattice_lines_config, nanotape_vacuum_config],
        direction=AxisEnum.y,
    )
    nanotape_builder = NanoTapeBuilder()
    nanotape_material = nanotape_builder.get_material(nanotape_config)
    vacuum_config = create_vacuum_config(
        size=vacuum_length,
        crystal=nanotape_material,
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
