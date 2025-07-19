from typing import List, Tuple, Optional

import numpy as np
from mat3ra.code.array_with_ids import ArrayWithIds
from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum

from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.interface import (
    CommensurateLatticeInterfaceAnalyzer,
    InterfaceAnalyzer,
    ZSLInterfaceAnalyzer,
)
from mat3ra.made.tools.analyze.interface.twisted_nanoribbons import TwistedNanoribbonsInterfaceAnalyzer
from .builders import (
    InterfaceBuilder,
)
from .configuration import (
    InterfaceConfiguration,
)
from .. import MaterialWithBuildMetadata
from ..slab.configurations import (
    SlabConfiguration,
    SlabStrainedSupercellConfiguration,
)
from ..vacuum.configuration import VacuumConfiguration
from ...calculate.calculators import InterfaceMaterialCalculator
from ...modify import interface_displace_part
from ...optimize import evaluate_calculator_on_xy_grid
from mat3ra.made.tools.analyze.slab import SlabMaterialAnalyzer


def create_simple_interface_between_slabs(
    substrate_slab: MaterialWithBuildMetadata,
    film_slab: MaterialWithBuildMetadata,
    gap: Optional[float] = None,
    vacuum: float = 10.0,
    xy_shift: List[float] = [0, 0],
) -> Material:
    """
    Create an interface between two slab materials with specified parameters.
    No strain matching is performed, the interface is created as a simple stack of two slabs.

    Args:
        substrate_slab (MaterialWithBuildMetadata): Substrate slab material.
        film_slab (MaterialWithBuildMetadata): Film slab material.
        gap (Optional[float]): Gap between the two materials, in Angstroms.
            Distance between top most atom of the substrate and bottom most atom of the film.
        vacuum (float): Size of the vacuum layer in Angstroms.
        xy_shift (List[float]): Shift in x and y directions, in Angstroms.

    Returns:
        Material: The interface material.
    """
    substrate_analyzer = SlabMaterialAnalyzer(material=substrate_slab)
    film_analyzer = SlabMaterialAnalyzer(material=film_slab)

    analyzer = InterfaceAnalyzer(
        substrate_slab_configuration=substrate_analyzer.build_configuration,
        film_slab_configuration=film_analyzer.build_configuration,
        substrate_build_parameters=substrate_analyzer.build_parameters,
        film_build_parameters=film_analyzer.build_parameters,
    )

    film_configuration = analyzer.film_strained_configuration
    substrate_configuration = analyzer.substrate_strained_configuration
    vacuum_configuration = VacuumConfiguration(
        size=vacuum,
    )

    config = InterfaceConfiguration(
        stack_components=[substrate_configuration, film_configuration, vacuum_configuration],
        xy_shift=xy_shift,
        gaps=ArrayWithIds.from_values([gap]),
    )
    builder = InterfaceBuilder()
    interface = builder.get_material(config)
    return interface


def create_zsl_interface(
    substrate_crystal: Material,
    film_crystal: Material,
    substrate_miller_indices: Tuple[int, int, int] = (0, 0, 1),
    film_miller_indices: Tuple[int, int, int] = (0, 0, 1),
    substrate_number_of_layers: int = 1,
    film_number_of_layers: int = 1,
    substrate_termination_formula: Optional[str] = None,
    film_termination_formula: Optional[str] = None,
    gap: Optional[float] = None,
    vacuum: float = 10.0,
    xy_shift: List[float] = [0, 0],
    max_area: float = 50.0,
    match_id: int = 0,
    max_area_ratio_tol: float = 0.09,
    max_length_tol: float = 0.03,
    max_angle_tol: float = 0.01,
) -> Material:
    """
    Create an interface between two materials using the ZSL algorithm to find matching supercells.

    Args:
        substrate_crystal (Material): Substrate bulk material.
        film_crystal (Material): Film bulk material.
        substrate_miller_indices (Tuple[int, int, int]): Miller indices for the substrate.
        film_miller_indices (Tuple[int, int, int]): Miller indices for the film.
        substrate_number_of_layers (int): Number of layers in the substrate.
        film_number_of_layers (int): Number of layers in the film.
        substrate_termination_formula (Optional[str]): Termination formula for the substrate.
        film_termination_formula (Optional[str]): Termination formula for the film.
        gap (Optional[float]): Gap between the two materials, in Angstroms.
            Distance between top most atom of the substrate and bottom most atom of the film.
        vacuum (float): Size of the vacuum layer in Angstroms.
        xy_shift (List[float]): Shift in x and y directions, in Angstroms.
        max_area (float): Maximum area for ZSL matching.
        match_id (int): ZSL match ID to use (0 for first match).
        max_area_ratio_tol (float): Area ratio tolerance for ZSL matching.
        max_length_tol (float): Length tolerance for ZSL matching.
        max_angle_tol (float): Angle tolerance for ZSL matching.

    Returns:
        Material: The ZSL interface material.
    """
    substrate_slab_config = SlabConfiguration.from_parameters(
        material_or_dict=substrate_crystal,
        miller_indices=substrate_miller_indices,
        number_of_layers=substrate_number_of_layers,
        vacuum=0,
        termination_formula=substrate_termination_formula,
    )
    film_slab_config = SlabConfiguration.from_parameters(
        material_or_dict=film_crystal,
        miller_indices=film_miller_indices,
        number_of_layers=film_number_of_layers,
        vacuum=0,
        termination_formula=film_termination_formula,
    )

    # Use ZSLInterfaceAnalyzer to get strained slab configurations
    analyzer = ZSLInterfaceAnalyzer(
        substrate_slab_configuration=substrate_slab_config,
        film_slab_configuration=film_slab_config,
        max_area=max_area,
        max_area_ratio_tol=max_area_ratio_tol,
        max_length_tol=max_length_tol,
        max_angle_tol=max_angle_tol,
    )

    selected_config = analyzer.get_strained_configuration_by_match_id(match_id=match_id)

    vacuum_configuration = VacuumConfiguration(
        size=vacuum,
    )
    interface_config = InterfaceConfiguration(
        stack_components=[
            selected_config.substrate_configuration,
            selected_config.film_configuration,
            vacuum_configuration,
        ],
        xy_shift=xy_shift,
        gaps=ArrayWithIds.from_values([gap]),
    )

    builder = InterfaceBuilder()
    return builder.get_material(interface_config)


def create_zsl_interface_between_slabs(
    substrate_slab: MaterialWithBuildMetadata,
    film_slab: MaterialWithBuildMetadata,
    gap: Optional[float] = None,
    vacuum: float = 10.0,
    xy_shift: Optional[List[float]] = None,
    max_area: float = 200.0,
    max_area_ratio_tol: float = 0.1,
    max_length_tol: float = 0.03,
    max_angle_tol: float = 0.01,
    match_id: int = 0,
) -> Material:
    """
    Create an interface between two slab materials using the ZSL algorithm.

    Args:
        substrate_slab (MaterialWithBuildMetadata): Substrate slab material.
        film_slab (MaterialWithBuildMetadata): Film slab material.
        gap (Optional[float]): Gap between the two materials, in Angstroms.
        vacuum (float): Size of the vacuum layer in Angstroms.
        xy_shift (Optional[List[float]]): Shift in x and y directions, in Angstroms. Defaults to [0, 0].
        max_area (float): Maximum area of the supercell.
        max_area_ratio_tol (float): Tolerance for the ratio of the areas of the two supercells.
        max_length_tol (float): Tolerance for the length of the supercell vectors.
        max_angle_tol (float): Tolerance for the angle between the supercell vectors.
        match_id (int): Index of the match to use from the list of found ZSL matches.

    Returns:
        Material: The interface material.
    """
    if xy_shift is None:
        xy_shift = [0, 0]

    substrate_slab_config = SlabMaterialAnalyzer(material=substrate_slab).build_configuration
    film_slab_config = SlabMaterialAnalyzer(material=film_slab).build_configuration

    analyzer = ZSLInterfaceAnalyzer(
        substrate_slab_configuration=substrate_slab_config,
        film_slab_configuration=film_slab_config,
        max_area=max_area,
        max_area_ratio_tol=max_area_ratio_tol,
        max_length_tol=max_length_tol,
        max_angle_tol=max_angle_tol,
    )
    interface_configurations = analyzer.get_strained_configurations()
    if not interface_configurations:
        raise ValueError("No ZSL match found for the given parameters.")

    if match_id >= len(interface_configurations):
        raise IndexError(f"match_id {match_id} is out of bounds for {len(interface_configurations)} matches found.")

    selected_config = interface_configurations[match_id]

    vacuum_configuration = VacuumConfiguration(
        size=vacuum,
    )

    interface_config = InterfaceConfiguration(
        stack_components=[
            selected_config.substrate_configuration,
            selected_config.film_configuration,
            vacuum_configuration,
        ],
        gaps=ArrayWithIds.from_values(values=[gap, gap] if gap is not None else []),
        xy_shift=xy_shift,
    )

    builder = InterfaceBuilder()
    return builder.get_material(interface_config)


def create_twisted_interface(
    material1: Material,
    material2: Material,
    angle: float = 0.0,
    vacuum_x: float = 5.0,
    vacuum_y: float = 5.0,
    gap: float = 3.0,
) -> Material:
    """
    Create a twisted interface between two nanoribbons.

    Args:
        material1 (Material): First nanoribbon material.
        material2 (Material): Second nanoribbon material.
        angle (float): Twist angle in degrees.
        vacuum_x (float): Vacuum along x on both sides, in Angstroms.
        vacuum_y (float): Vacuum along y on both sides, in Angstroms.
        gap (float): Gap between the nanoribbons in Angstroms.

    Returns:
        Material: The twisted interface material.
    """
    slab1 = SlabConfiguration.from_parameters(
        material_or_dict=material1,
        miller_indices=(0, 0, 1),
        number_of_layers=1,
        vacuum=0.0,
    )
    slab2 = SlabConfiguration.from_parameters(
        material_or_dict=material2,
        miller_indices=(0, 0, 1),
        number_of_layers=1,
        vacuum=0.0,
    )
    analyzer = TwistedNanoribbonsInterfaceAnalyzer(
        substrate_slab_configuration=slab1,
        film_slab_configuration=slab2,
        angle=angle,
        vacuum_x=vacuum_x,
        vacuum_y=vacuum_y,
    )
    processed_slab1 = analyzer.substrate_nanoribbon_configuration
    processed_slab2 = analyzer.film_nanoribbon_configuration

    configuration = InterfaceConfiguration(
        stack_components=[processed_slab1, processed_slab2],
        gaps=ArrayWithIds.from_values([gap, gap]),
        direction=AxisEnum.z,
    )
    builder = InterfaceBuilder()
    return builder.get_material(configuration)


def get_commensurate_strained_configurations(
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
) -> Tuple[List[SlabStrainedSupercellConfiguration], float]:
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

    Returns:
        Tuple[List[SlabStrainedSupercellConfiguration], float]:
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
    actual_angle = match_holders[match_id].actual_angle
    return [selected_config.substrate_configuration, selected_config.film_configuration], actual_angle


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
    strained_configs, _ = get_commensurate_strained_configurations(
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
    )

    interface_config = InterfaceConfiguration(
        stack_components=strained_configs,
        gaps=ArrayWithIds.from_values([gap, gap]),
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
