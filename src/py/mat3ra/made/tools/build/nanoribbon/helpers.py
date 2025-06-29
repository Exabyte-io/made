from typing import Tuple, Optional, Union

from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.nanoribbon_analyzer import NanoribbonAnalyzer
from mat3ra.made.tools.analyze.nanotape_analyzer import NanoTapeAnalyzer
from mat3ra.made.tools.build.slab.entities import Termination
from . import EdgeTypes
from .builders import NanoribbonBuilder, NanoTapeBuilder, NanoTapeBuilderParameters, NanoribbonBuilderParameters
from .configuration import (
    get_miller_indices_from_edge_type,
)


def create_nanotape(
    material: Material,
    miller_indices_uv: Optional[Tuple[int, int]] = None,
    edge_type: Optional[EdgeTypes] = EdgeTypes.zigzag,
    width: int = 2,
    length: int = 2,
    vacuum_width: float = 10.0,
    use_rectangular_cell: bool = True,
) -> Material:
    """
    Create a nanotape from a monolayer material.

    Args:
        material: The monolayer material to create the nanotape from.
        miller_indices_uv: The (u,v) Miller indices.
        edge_type: Edge type string ("zigzag"/"armchair"). Optional if miller_indices_uv is provided.
        width: The width of the nanotape in number of unit cells.
        length: The length of the nanotape in number of unit cells.
        vacuum_width: The width of the vacuum region in Angstroms (cartesian).
        use_rectangular_cell: Whether the nanotape is rectangular.

    Returns:
        The nanotape material.
    """
    if miller_indices_uv is None and edge_type is None:
        raise ValueError("Either miller_indices_uv or edge_type must be provided")

    if miller_indices_uv is None and edge_type is not None:
        miller_indices_uv = get_miller_indices_from_edge_type(edge_type)

    analyzer = NanoTapeAnalyzer(material=material, miller_indices_uv=miller_indices_uv)
    configuration = analyzer.get_configuration(
        width=width,
        length=length,
        vacuum_width=vacuum_width,
    )

    build_parameters = NanoTapeBuilderParameters(use_rectangular_lattice=use_rectangular_cell)
    builder = NanoTapeBuilder(build_parameters=build_parameters)
    return builder.get_material(configuration)


def create_nanoribbon(
    material: Material,
    miller_indices_uv: Optional[Union[Tuple[int, int], str]] = None,
    edge_type: Optional[EdgeTypes] = EdgeTypes.zigzag,
    width: int = 2,
    length: int = 2,
    vacuum_width: float = 10.0,
    vacuum_length: float = 0.0,
    use_rectangular_cell: bool = True,
    termination: Optional[Termination] = None,
) -> Material:
    """
    Creates a nanoribbon material from a monolayer material.

    Args:
        material: The monolayer material to create the nanoribbon from (assumes vacuum is present).
        miller_indices_uv: The (u,v) Miller indices for the nanoribbon direction.
        edge_type: Edge type string ("zigzag"/"armchair"). Optional if miller_indices_uv is provided.
        width: The width of the nanoribbon in number of unit cells.
        length: The length of the nanoribbon in number of unit cells.
        vacuum_width: The width of the vacuum region in Angstroms (cartesian).
        vacuum_length: The length of the vacuum region in Angstroms (cartesian).
        use_rectangular_cell: Whether the nanoribbon is rectangular.
        termination: The termination to use for the nanoribbon. If None, uses default termination.

    Returns:
        Material: The generated nanoribbon material.
    """
    if miller_indices_uv is None and edge_type is None:
        raise ValueError("Either miller_indices_uv or edge_type must be provided")

    if miller_indices_uv is None and edge_type is not None:
        miller_indices_uv = get_miller_indices_from_edge_type(edge_type)

    analyzer = NanoribbonAnalyzer(material=material, miller_indices_uv=miller_indices_uv)
    configuration = analyzer.get_configuration(
        width=width,
        length=length,
        vacuum_width=vacuum_width,
        vacuum_length=vacuum_length,
    )

    build_parameters = NanoribbonBuilderParameters(use_rectangular_lattice=use_rectangular_cell)
    builder = NanoribbonBuilder(build_parameters=build_parameters)
    nanoribbon = builder.get_material(configuration)

    return nanoribbon
