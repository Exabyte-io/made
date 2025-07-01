from typing import Optional, Union, List, Tuple

from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.interface.grain_boundary import GrainBoundaryAnalyzer
from .builders import GrainBoundaryBuilder
from .configuration import GrainBoundaryConfiguration
from mat3ra.made.tools.build.vacuum.configuration import VacuumConfiguration
from mat3ra.made.tools.build.slab.configurations import SlabStrainedSupercellWithGapConfiguration


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

    # Create grain boundary configuration using the new structure
    # Convert translation_vector to xy_shift (y,z components for x-direction stacking)
    xy_shift = [translation_vector[1], translation_vector[2]]
    
    # Wrap configurations with gap if specified
    if gap > 0:
        substrate_config = SlabStrainedSupercellWithGapConfiguration(
            **strained_config.substrate_configuration.to_dict(),
            gap=gap
        )
        film_config = SlabStrainedSupercellWithGapConfiguration(
            **strained_config.film_configuration.to_dict(),
            gap=gap
        )
    else:
        substrate_config = strained_config.substrate_configuration
        film_config = strained_config.film_configuration
    
    gb_config = GrainBoundaryConfiguration(
        stack_components=[substrate_config, film_config],
        direction=AxisEnum.x,
        xy_shift=xy_shift,
    )

    # Build the grain boundary
    builder = GrainBoundaryBuilder()
    return builder.get_material(gb_config)


def create_grain_boundary(
    configuration: Union[GrainBoundaryConfiguration],
    builder: Optional[Union[GrainBoundaryBuilder]] = None,
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
        builder = GrainBoundaryBuilder()

    return builder.get_material(configuration)
