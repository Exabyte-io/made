from typing import Optional, List, Tuple

from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum

from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.interface.commensurate import CommensurateLatticeInterfaceAnalyzer
from .builders import GrainBoundaryLinearBuilder, SlabGrainBoundaryBuilder, SlabGrainBoundaryBuilderParameters
from .configuration import GrainBoundaryLinearConfiguration, SlabGrainBoundaryConfiguration
from ..slab.configurations import SlabConfiguration, SlabStrainedSupercellWithGapConfiguration


def create_grain_boundary_planar(
    phase_1_material: Material,
    phase_2_material: Optional[Material] = None,
    phase_1_miller_indices: Tuple[int, int, int] = (0, 0, 1),
    phase_2_miller_indices: Tuple[int, int, int] = (0, 0, 1),
    phase_1_thickness: int = 1,
    phase_2_thickness: int = 1,
    gap: float = 3.0,
    slab_vacuum: float = 1.0,
    slab_xy_supercell_matrix: List[List[int]] = [[1, 0], [0, 1]],
) -> Material:
    """
    Create a planar grain boundary between two materials with different orientations.

    This function creates a grain boundary by:
    1. Creating slab configurations for both phases
    2. Creating an interface between the phases
    3. Rotating the interface by 90 degrees to create a planar grain boundary

    Args:
        phase_1_material (Material): First phase material.
        phase_2_material (Optional[Material]): Second phase material. If None, uses phase_1_material.
        phase_1_miller_indices (Tuple[int, int, int]): Miller indices for phase 1.
        phase_2_miller_indices (Tuple[int, int, int]): Miller indices for phase 2.
        phase_1_thickness (int): Number of layers for phase 1.
        phase_2_thickness (int): Number of layers for phase 2.
        gap (float): Gap between phases in Angstroms.
        slab_vacuum (float): Vacuum for the final slab in Angstroms.
        slab_xy_supercell_matrix (List[List[int]]): Supercell matrix for the final slab.

    Returns:
        Material: The planar grain boundary material.
    """
    if phase_2_material is None:
        phase_2_material = phase_1_material

    # Create slab configurations for both phases
    phase_1_configuration = SlabConfiguration.from_parameters(
        material_or_dict=phase_1_material,
        miller_indices=phase_1_miller_indices,
        number_of_layers=phase_1_thickness,
        vacuum=0.0,
    )
    phase_2_configuration = SlabConfiguration.from_parameters(
        material_or_dict=phase_2_material,
        miller_indices=phase_2_miller_indices,
        number_of_layers=phase_2_thickness,
        vacuum=0.0,
    )

    # Get terminations
    termination1 = phase_1_configuration.get_terminations()[0]
    termination2 = phase_2_configuration.get_terminations()[0]

    # Create slab configuration for the final grain boundary
    slab_config = SlabConfiguration.from_parameters(
        material_or_dict=phase_1_material,
        miller_indices=phase_1_miller_indices,
        number_of_layers=phase_1_thickness,
        vacuum=slab_vacuum,
        xy_supercell_matrix=slab_xy_supercell_matrix,
    )

    # Create grain boundary configuration
    config = SlabGrainBoundaryConfiguration(
        phase_1_configuration=phase_1_configuration,
        phase_2_configuration=phase_2_configuration,
        phase_1_termination=termination1,
        phase_2_termination=termination2,
        gap=gap,
        slab_configuration=slab_config,
    )

    # Build the grain boundary
    builder_params = SlabGrainBoundaryBuilderParameters()
    builder = SlabGrainBoundaryBuilder(build_parameters=builder_params)
    grain_boundary = builder.get_material(config)

    return grain_boundary


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
    3. Creating strained configurations for both phases
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
    # Create slab configuration from the material
    slab_config = SlabConfiguration.from_parameters(
        material_or_dict=material,
        miller_indices=miller_indices,
        number_of_layers=number_of_layers,
        vacuum=vacuum,
    )

    # Find commensurate lattice matches using the analyzer
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

    # Get the first match
    selected_config = analyzer.get_strained_configuration_by_match_id(0)

    # Create strained configurations with gap
    substrate_config_with_gap = SlabStrainedSupercellWithGapConfiguration(
        **selected_config.substrate_configuration.to_dict(), gap=gap
    )
    film_config_with_gap = SlabStrainedSupercellWithGapConfiguration(
        **selected_config.film_configuration.to_dict(), gap=gap
    )

    # Create grain boundary configuration
    grain_boundary_config = GrainBoundaryLinearConfiguration(
        stack_components=[substrate_config_with_gap, film_config_with_gap],
        direction=direction,
        gap=gap,
    )

    # Build the grain boundary using default parameters
    builder = GrainBoundaryLinearBuilder()
    grain_boundary = builder.get_material(grain_boundary_config)

    return grain_boundary
