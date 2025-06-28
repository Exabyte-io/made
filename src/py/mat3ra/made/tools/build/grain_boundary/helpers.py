from typing import Optional, Union, List, Tuple

from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.interface.grain_boundary import GrainBoundaryAnalyzer
from mat3ra.made.tools.analyze.interface.utils.holders import MatchedSubstrateFilmConfigurationHolder

from .builders import GrainBoundaryBuilder, GrainBoundaryWithVacuumBuilder
from .configuration import GrainBoundaryConfiguration, GrainBoundaryWithVacuumConfiguration


def create_grain_boundary_planar(
    phase_1_material: Material,
    phase_2_material: Material,
    phase_1_miller_indices: Tuple[int, int, int],
    phase_2_miller_indices: Tuple[int, int, int],
    phase_1_thickness: int = 1,
    phase_2_thickness: int = 1,
    translation_vector: List[float] = [0.0, 0.0, 0.0],
    gap: float = 3.0,
    match_id: int = 0,
    max_area: float = 50.0,
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
    # Create analyzer to get configurations
    analyzer = GrainBoundaryAnalyzer(
        phase_1_material=phase_1_material,
        phase_2_material=phase_2_material,
        phase_1_miller_indices=phase_1_miller_indices,
        phase_2_miller_indices=phase_2_miller_indices,
        phase_1_thickness=phase_1_thickness,
        phase_2_thickness=phase_2_thickness,
        max_area=max_area,
        max_area_ratio_tol=max_area_ratio_tol,
        max_length_tol=max_length_tol,
        max_angle_tol=max_angle_tol,
    )

    # Get the strained configuration for the specified match
    strained_config = analyzer.get_grain_boundary_configuration_by_match_id(match_id)

    # Create grain boundary configuration
    gb_config = GrainBoundaryConfiguration(
        phase_1_configuration=strained_config,
        phase_2_configuration=strained_config,
        translation_vector=translation_vector,
        gap=gap,
    )

    # Build the grain boundary
    builder = GrainBoundaryBuilder()
    return builder.get_material(gb_config)


def create_grain_boundary_planar_with_vacuum(
    phase_1_material: Material,
    phase_2_material: Material,
    phase_1_miller_indices: Tuple[int, int, int],
    phase_2_miller_indices: Tuple[int, int, int],
    slab_miller_indices: Tuple[int, int, int],
    phase_1_thickness: int = 1,
    phase_2_thickness: int = 1,
    slab_thickness: int = 1,
    translation_vector: List[float] = [0.0, 0.0, 0.0],
    gap: float = 3.0,
    vacuum: float = 10.0,
    termination_formula: Optional[str] = None,
    xy_supercell_matrix: Optional[List[List[int]]] = None,
    match_id: int = 0,
    max_area: float = 50.0,
    max_area_ratio_tol: float = 0.09,
    max_length_tol: float = 0.03,
    max_angle_tol: float = 0.01,
) -> Material:
    """
    Create a planar grain boundary with vacuum (slab).

    Args:
        phase_1_material: First phase material
        phase_2_material: Second phase material
        phase_1_miller_indices: Miller indices for phase 1
        phase_2_miller_indices: Miller indices for phase 2
        slab_miller_indices: Miller indices for the slab surface
        phase_1_thickness: Number of layers for phase 1
        phase_2_thickness: Number of layers for phase 2
        slab_thickness: Number of layers for the slab
        translation_vector: Relative shift between phases [x, y, z]
        gap: Gap between phases in Angstroms
        vacuum: Vacuum size in Angstroms
        termination_formula: Termination formula for the slab
        xy_supercell_matrix: Supercell matrix for xy plane
        match_id: ZSL match ID to use (0 for first match)
        max_area: Maximum area for ZSL matching
        max_area_ratio_tol: Area ratio tolerance for ZSL matching
        max_length_tol: Length tolerance for ZSL matching
        max_angle_tol: Angle tolerance for ZSL matching

    Returns:
        Material: The grain boundary slab material
    """
    # First create the planar grain boundary
    grain_boundary_material = create_grain_boundary_planar(
        phase_1_material=phase_1_material,
        phase_2_material=phase_2_material,
        phase_1_miller_indices=phase_1_miller_indices,
        phase_2_miller_indices=phase_2_miller_indices,
        phase_1_thickness=phase_1_thickness,
        phase_2_thickness=phase_2_thickness,
        translation_vector=translation_vector,
        gap=gap,
        match_id=match_id,
        max_area=max_area,
        max_area_ratio_tol=max_area_ratio_tol,
        max_length_tol=max_length_tol,
        max_angle_tol=max_angle_tol,
    )

    # Create configuration for grain boundary with vacuum
    gb_vacuum_config = GrainBoundaryWithVacuumConfiguration.from_grain_boundary_material(
        grain_boundary_material=grain_boundary_material,
        miller_indices=slab_miller_indices,
        number_of_layers=slab_thickness,
        vacuum=vacuum,
        termination_formula=termination_formula,
        xy_supercell_matrix=xy_supercell_matrix,
    )

    # Build the grain boundary with vacuum
    builder = GrainBoundaryWithVacuumBuilder()
    return builder.get_material(gb_vacuum_config)


def create_grain_boundary(
    configuration: Union[GrainBoundaryConfiguration, GrainBoundaryWithVacuumConfiguration],
    builder: Optional[Union[GrainBoundaryBuilder, GrainBoundaryWithVacuumBuilder]] = None,
) -> Material:
    """
    Create a grain boundary according to provided configuration with selected builder.

    Args:
        configuration: The configuration of the grain boundary
        builder: The builder to use for creating the grain boundary

    Returns:
        Material: The material with the grain boundary
    """
    if builder is None:
        if isinstance(configuration, GrainBoundaryConfiguration):
            builder = GrainBoundaryBuilder()
        else:
            builder = GrainBoundaryWithVacuumBuilder()

    return builder.get_material(configuration)
