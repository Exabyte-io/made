from typing import Tuple, Optional, Type, TypeVar

from mat3ra.made.material import Material
from mat3ra.made.tools.build import BaseBuilder, BaseBuilderParameters
from mat3ra.made.tools.build.slab.entities import Termination
from . import NanoribbonConfiguration
from .builders import NanoribbonBuilder, NanoribbonBuilderParameters
from ..lattice_lines.configuration import EdgeTypes, get_miller_indices_from_edge_type
from ..nanotape import NanoTapeConfiguration
from ..nanotape.builders import NanoTapeBuilder, NanoTapeBuilderParameters

T = TypeVar("T", bound=BaseBuilder)
P = TypeVar("P", bound=BaseBuilderParameters)


def _resolve_miller_indices(
    miller_indices_uv: Optional[Tuple[int, int]], edge_type: Optional[EdgeTypes]
) -> Tuple[int, int]:
    """Resolve Miller indices from either direct input or edge type."""
    if miller_indices_uv is None and edge_type is None:
        raise ValueError("Either miller_indices_uv or edge_type must be provided")

    if miller_indices_uv is None and edge_type is not None:
        miller_indices_uv = get_miller_indices_from_edge_type(edge_type)

    assert miller_indices_uv is not None  # Type checker hint
    return miller_indices_uv


def _create_material_with_analyzer_and_builder(
    material: Material,
    miller_indices_uv: Tuple[int, int],
    analyzer_class: Type,
    builder_class: Type[T],
    build_parameters_class: Type[P],
    use_rectangular_cell: bool,
    **analyzer_kwargs,
) -> Material:
    """Generic function to create material using analyzer and builder pattern."""
    analyzer = analyzer_class(material=material, miller_indices_uv=miller_indices_uv)
    configuration = analyzer.get_configuration(**analyzer_kwargs)

    build_parameters = build_parameters_class(use_rectangular_lattice=use_rectangular_cell)
    builder = builder_class(build_parameters=build_parameters)
    return builder.get_material(configuration)


def create_nanoribbon(
    material: Material,
    miller_indices_uv: Optional[Tuple[int, int]] = None,
    edge_type: Optional[EdgeTypes] = EdgeTypes.zigzag,
    width: int = 2,
    length: int = 2,
    vacuum_width: float = 10.0,
    vacuum_length: float = 10.0,
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
    uv = _resolve_miller_indices(miller_indices_uv, edge_type)

    is_tape = vacuum_length is None
    if is_tape:
        builder_cls = NanoTapeBuilder
        params_cls = NanoTapeBuilderParameters
        config = NanoTapeConfiguration.from_parameters(
            material=material,
            miller_indices_uv=uv,
            width=width,
            length=length,
            vacuum_width=vacuum_width,
            termination=termination,
        )
    else:
        builder_cls = NanoribbonBuilder
        params_cls = NanoribbonBuilderParameters
        config = NanoribbonConfiguration.from_parameters(
            material=material,
            miller_indices_uv=uv,
            width=width,
            length=length,
            vacuum_width=vacuum_width,
            vacuum_length=vacuum_length,
            termination=termination,
        )

    builder = builder_cls(build_parameters=params_cls(use_rectangular_cell=use_rectangular_cell))
    return builder.get_material(config)
