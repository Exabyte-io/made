from typing import List, Tuple, Optional

import numpy as np

from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.interface import CommensurateLatticeInterfaceAnalyzer
from .builders import (
    InterfaceBuilder,
    InterfaceBuilderParameters,
)
from .configuration import (
    InterfaceConfiguration,
)
from ..slab.configurations import (
    SlabConfiguration,
    SlabStrainedSupercellWithGapConfiguration,
)
from ...calculate.calculators import InterfaceMaterialCalculator
from ...modify import interface_displace_part
from ...optimize import evaluate_calculator_on_xy_grid


def get_commensurate_strained_configs(
    material: Material,
    target_angle: float,
    angle_tolerance: float,
    max_repetition_int: Optional[int],
    limit_max_int: int,
    return_first_match: bool,
    miller_indices: Tuple[int, int, int],
    number_of_layers: int,
    vacuum: float,
    match_id: int = 0,
    direction: AxisEnum = AxisEnum.z,
    gap: float = 3.0,
) -> Tuple[List[SlabStrainedSupercellWithGapConfiguration], float]:
    """
    Get strained configurations for commensurate lattice matching.

    This function creates strained configurations for both substrate and film phases
    based on commensurate lattice matching at a specified twist angle.

    Args:
        material (Material): The material to create configurations from.
        target_angle (float): The target twist angle in degrees.
        angle_tolerance (float): Tolerance for matching angles in degrees.
        max_repetition_int (Optional[int]): Maximum integer for supercell matrix elements.
        limit_max_int (int): Limit for maximum integer to search for supercell matrices.
        return_first_match (bool): Whether to return the first match or all matches.
        miller_indices (Tuple[int, int, int]): Miller indices for the slab surface.
        number_of_layers (int): Number of atomic layers in the slab.
        vacuum (float): Size of the vacuum layer in Angstroms.
        match_id (int): ID of the match to use (0 for first match).
        direction (AxisEnum): Direction for the gap configuration.
        gap (float): The gap between the two phases in Angstroms.

    Returns:
        Tuple[List[SlabStrainedSupercellWithGapConfiguration], float]:
            List of strained configurations [substrate, film] and the actual angle.

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

    selected_config = analyzer.get_strained_configuration_by_match_id(match_id)

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

    actual_angle = match_holders[match_id].actual_angle
    return [substrate_config, film_config], actual_angle


def create_interface(
    configuration: InterfaceConfiguration,
) -> Material:
    builder = InterfaceBuilder()
    return builder.get_material(configuration)


def create_commensurate_interface(
    material: Material,
    target_angle: float = 0.0,
    angle_tolerance: float = 0.1,
    max_repetition_int: Optional[int] = None,
    limit_max_int: int = 20,
    return_first_match: bool = True,
    direction: AxisEnum = AxisEnum.z,
    gap: float = 3.0,
    miller_indices: Tuple[int, int, int] = (0, 0, 1),
    number_of_layers: int = 1,
    vacuum: float = 0.0,
    match_id: int = 0,
) -> Material:
    """
    Create a commensurate lattice interface (twisted interface) from a material with specified twist angle.

    This function creates a commensurate interface by:
    1. Creating a slab configuration from the material
    2. Finding commensurate lattice matches at the target angle using CommensurateLatticeInterfaceAnalyzer
    3. Creating strained configurations with directional gaps for both phases
    4. Stacking them along the specified direction

    Args:
        material (Material): The material to create the interface from.
        target_angle (float): The target twist angle in degrees.
        angle_tolerance (float): Tolerance for matching angles in degrees.
        max_repetition_int (Optional[int]): Maximum integer for supercell matrix elements.
        limit_max_int (int): Limit for maximum integer to search for supercell matrices.
        return_first_match (bool): Whether to return the first match or all matches.
        direction (AxisEnum): Direction along which to stack components (x, y, or z).
        gap (float): The gap between the two phases in Angstroms.
        miller_indices (Tuple[int, int, int]): Miller indices for the slab surface.
        number_of_layers (int): Number of atomic layers in the slab.
        vacuum (float): Size of the vacuum layer in Angstroms.
        match_id (int): ID of the match to use (0 for first match).

    Returns:
        Material: The commensurate interface material.

    Raises:
        ValueError: If no commensurate lattice matches are found.
    """
    strained_configs, _ = get_commensurate_strained_configs(
        material=material,
        target_angle=target_angle,
        angle_tolerance=angle_tolerance,
        max_repetition_int=max_repetition_int,
        limit_max_int=limit_max_int,
        return_first_match=return_first_match,
        miller_indices=miller_indices,
        number_of_layers=number_of_layers,
        vacuum=vacuum,
        match_id=match_id,
        direction=direction,
        gap=gap,
    )

    interface_config = InterfaceConfiguration(
        stack_components=strained_configs,
        direction=direction,
    )

    builder = InterfaceBuilder()
    return builder.get_material(interface_config)


def get_optimal_film_displacement(
    material: Material,
    grid_size_xy: Tuple[int, int] = (10, 10),
    grid_offset_position: List[float] = [0, 0],
    grid_range_x=(-0.5, 0.5),
    grid_range_y=(-0.5, 0.5),
    use_cartesian_coordinates=False,
    calculator: InterfaceMaterialCalculator = InterfaceMaterialCalculator(),
):
    """
    Calculate the optimal displacement in of the film to minimize the interaction energy
        between the film and the substrate. The displacement is calculated on a grid.

    This function evaluates the interaction energy between the film and substrate
    over a specified grid of (x,y) displacements. It returns the displacement vector that
    results in the minimum interaction energy.

    Args:
        material (Material): The interface Material object.
        grid_size_xy (Tuple[int, int]): The size of the grid to search for the optimal displacement.
        grid_offset_position (List[float]): The offset position of the grid.
        grid_range_x (Tuple[float, float]): The range of the grid in x.
        grid_range_y (Tuple[float, float]): The range of the grid in y.
        use_cartesian_coordinates (bool): Whether to use Cartesian coordinates.
        calculator (InterfaceMaterialCalculator): The calculator to use for the calculation of the interaction energy.

    Returns:
        List[float]: The optimal displacement vector.

    """
    xy_matrix, results_matrix = evaluate_calculator_on_xy_grid(
        material=material,
        calculator_function=calculator.get_energy,
        modifier=interface_displace_part,
        modifier_parameters={},
        grid_size_xy=grid_size_xy,
        grid_offset_position=grid_offset_position,
        grid_range_x=grid_range_x,
        grid_range_y=grid_range_y,
        use_cartesian_coordinates=use_cartesian_coordinates,
    )
    min_index = np.unravel_index(np.argmin(results_matrix), results_matrix.shape)

    optimal_x = xy_matrix[0][min_index[0]]
    optimal_y = xy_matrix[1][min_index[1]]

    return [optimal_x, optimal_y, 0]
