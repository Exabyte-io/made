from typing import Optional, List, Tuple

from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum

from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.interface.commensurate import CommensurateLatticeInterfaceAnalyzer
from mat3ra.made.tools.analyze.interface.grain_boundary import GrainBoundaryAnalyzer
from .builders import GrainBoundaryBuilder
from .builders import GrainBoundaryLinearBuilder
from .configuration import GrainBoundaryConfiguration
from .configuration import GrainBoundaryLinearConfiguration
from ..slab.configurations import (
    SlabConfiguration,
    SlabStrainedSupercellWithGapConfiguration,
)


def create_grain_boundary_planar(
    phase_1_material: Material,
    phase_2_material: Optional[Material] = None,
    phase_1_miller_indices: Tuple[int, int, int] = (0, 0, 1),
    phase_2_miller_indices: Tuple[int, int, int] = (0, 0, 1),
    phase_1_thickness: int = 1,
    phase_2_thickness: int = 1,
    translation_vector: List[float] = [0.0, 0.0],
    gap: float = 3.0,
    max_area: float = 50.0,
    match_id: int = 0,
    max_area_ratio_tol: float = 0.09,
    max_length_tol: float = 0.03,
    max_angle_tol: float = 0.01,
) -> Material:
    """
    Create a planar grain boundary between two materials with different orientations.

    Args:
        phase_1_material: First phase material
        phase_2_material: Second phase material
        phase_1_miller_indices: Miller indices for phase 1
        phase_2_miller_indices: Miller indices for phase 2
        phase_1_thickness: Number of layers for phase 1
        phase_2_thickness: Number of layers for phase 2
        translation_vector: Relative shift between phases [x, y, z]
        gap: Gap between phases in Angstroms
        match_id: ZSL match ID to use (0 for first match)
        max_area: Maximum area for ZSL matching
        max_area_ratio_tol: Area ratio tolerance for ZSL matching
        max_length_tol: Length tolerance for ZSL matching
        max_angle_tol: Angle tolerance for ZSL matching

    Returns:
        Material: The grain boundary material
    """
    analyzer = GrainBoundaryAnalyzer(
        phase_1_material=phase_1_material,
        phase_2_material=phase_2_material if phase_2_material else phase_1_material,
        phase_1_miller_indices=phase_1_miller_indices,
        phase_2_miller_indices=phase_2_miller_indices,
        phase_1_thickness=phase_1_thickness,
        phase_2_thickness=phase_2_thickness,
        max_area=max_area,
        max_area_ratio_tol=max_area_ratio_tol,
        max_length_tol=max_length_tol,
        max_angle_tol=max_angle_tol,
    )

    strained_config = analyzer.get_grain_boundary_configuration_by_match_id(match_id)

    gb_config = GrainBoundaryConfiguration.from_parameters(
        phase_1_configuration=strained_config.substrate_configuration,
        phase_2_configuration=strained_config.film_configuration,
        xy_shift=translation_vector,
        gap=gap,
    )

    builder = GrainBoundaryBuilder()
    return builder.get_material(gb_config)


def create_grain_boundary_linear(
    material: Material,
    target_angle: float = 0.0,
    angle_tolerance: float = 0.1,
    max_repetition_int: Optional[int] = None,
    limit_max_int: int = 20,
    return_first_match: bool = True,
    direction: AxisEnum = AxisEnum.x,
    gap: float = 3.0,
    miller_indices: Tuple[int, int, int] = (0, 0, 1),
    number_of_layers: int = 1,
    vacuum: float = 0.0,
) -> Material:
    """
    Create a linear grain boundary from a material with specified twist angle.

    This function creates a grain boundary by:
    1. Creating a slab configuration from the material
    2. Finding commensurate lattice matches at the target angle using CommensurateLatticeInterfaceAnalyzer
    3. Creating strained configurations with directional gaps for both phases
    4. Stacking them along the specified direction

    Args:
        material (Material): The material to create the grain boundary from.
        target_angle (float): The target twist angle in degrees.
        angle_tolerance (float): Tolerance for matching angles in degrees.
        max_repetition_int (Optional[int]): Maximum integer for supercell matrix elements.
        limit_max_int (int): Limit for maximum integer to search for supercell matrices.
        return_first_match (bool): Whether to return the first match or all matches.
        direction (AxisEnum): Direction along which to stack components (x or y).
        gap (float): The gap between the two phases in Angstroms.
        miller_indices (Tuple[int, int, int]): Miller indices for the slab surface.
        number_of_layers (int): Number of atomic layers in the slab.
        vacuum (float): Size of the vacuum layer in Angstroms.

    Returns:
        Material: The grain boundary material.

    Raises:
        ValueError: If no commensurate lattice matches are found.
    """
    slab_config = SlabConfiguration.from_parameters(
        material_or_dict=material,
        miller_indices=miller_indices,
        number_of_layers=number_of_layers,
        vacuum=vacuum,
    )

    analyzer = CommensurateLatticeInterfaceAnalyzer(
        substrate_slab_configuration=slab_config,
        target_angle=target_angle,
        angle_tolerance=angle_tolerance,
        max_supercell_matrix_int=max_repetition_int,
        limit_max_int=limit_max_int,
        return_first_match=return_first_match,
    )

    match_holders = analyzer.commensurate_lattice_match_holders
    if not match_holders:
        raise ValueError(f"No commensurate lattice matches found for angle {target_angle}°")

    selected_config = analyzer.get_strained_configuration_by_match_id(0)

    substrate_config = SlabStrainedSupercellWithGapConfiguration(
        **selected_config.substrate_configuration.to_dict(),
        gap=gap,
        gap_direction=direction,
    )
    film_config = SlabStrainedSupercellWithGapConfiguration(
        **selected_config.film_configuration.to_dict(),
        gap=gap,
        gap_direction=direction,
    )

    actual_angle = getattr(match_holders[0], "angle", target_angle)
    grain_boundary_config = GrainBoundaryLinearConfiguration(
        stack_components=[substrate_config, film_config],
        direction=direction,
        gap=gap,
        actual_angle=actual_angle,
    )

    builder = GrainBoundaryLinearBuilder()
    grain_boundary = builder.get_material(grain_boundary_config)

    return grain_boundary
