from typing import Tuple, Optional

from mat3ra.made.material import Material
from mat3ra.made.tools.build.slab.entities import Termination
from . import NanoTapeConfiguration
from .builders import NanoTapeBuilder, NanoTapeBuilderParameters
from ..lattice_lines.configuration import EdgeTypes, get_miller_indices_from_edge_type


def create_nanotape(
    material: Material,
    miller_indices_2d: Optional[Tuple[int, int]] = None,
    edge_type: EdgeTypes = EdgeTypes.zigzag,
    width: int = 2,
    length: int = 2,
    vacuum_width: float = 10.0,
    use_rectangular_lattice: bool = True,
    termination: Optional[Termination] = None,
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
        termination: The termination to use for the nanotape. If None, uses default termination.

    Returns:
        Material: The generated nanotape material.
    """
    if miller_indices_2d is None:
        miller_indices_2d = get_miller_indices_from_edge_type(edge_type)

    config = NanoTapeConfiguration.from_parameters(
        material=material,
        miller_indices_2d=miller_indices_2d,
        width=width,
        length=length,
        vacuum_width=vacuum_width,
        termination=termination,
    )

    builder = NanoTapeBuilder(
        build_parameters=NanoTapeBuilderParameters(use_rectangular_lattice=use_rectangular_lattice)
    )
    return builder.get_material(config)
