from typing import Union, List, Optional, Tuple

import numpy as np

from mat3ra.made.material import Material
from ...calculate.calculators import InterfaceMaterialCalculator
from ...modify import interface_displace_part
from ...optimize import evaluate_calculator_on_xy_grid
from .builders import (
    SimpleInterfaceBuilder,
    SimpleInterfaceBuilderParameters,
    ZSLStrainMatchingParameters,
    ZSLStrainMatchingInterfaceBuilder,
    ZSLStrainMatchingInterfaceBuilderParameters,
)
from .configuration import InterfaceConfiguration


def create_interfaces(
    builder: Union[SimpleInterfaceBuilder, ZSLStrainMatchingInterfaceBuilder], configuration: InterfaceConfiguration
) -> List[Material]:
    return builder.get_materials(configuration)


def create_interface(
    configuration: InterfaceConfiguration,
    builder: Optional[Union[SimpleInterfaceBuilder, ZSLStrainMatchingInterfaceBuilder]] = None,
) -> Material:
    if builder is None:
        builder = SimpleInterfaceBuilder(build_parameters=SimpleInterfaceBuilderParameters())
    return builder.get_material(configuration)


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
