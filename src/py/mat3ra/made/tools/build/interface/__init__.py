from typing import Union, List, Optional, Tuple

import numpy as np

from mat3ra.made.material import Material
from ...calculate import calculate_norm_of_distances
from ...modify import displace_interface_part
from ...analyze import calculate_on_xy_grid
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
):
    results_matrix = calculate_on_xy_grid(
        material,
        modifier=displace_interface_part,
        modifier_parameters={},
        calculator=calculate_norm_of_distances,
        calculator_parameters={"shadowing_radius": 2.5},
        grid_size_xy=grid_size_xy,
        grid_offset_position=grid_offset_position,
        grid_range_x=grid_range_x,
        grid_range_y=grid_range_y,
        use_cartesian_coordinates=use_cartesian_coordinates,
    )
    min_index = np.unravel_index(np.argmin(results_matrix), results_matrix.shape)
    min_value = results_matrix[min_index]

    # Calculate the corresponding x, y coordinates
    x_values = np.linspace(grid_range_x[0], grid_range_x[1], grid_size_xy[0]) + grid_offset_position[0]
    y_values = np.linspace(grid_range_y[0], grid_range_y[1], grid_size_xy[1]) + grid_offset_position[1]

    optimal_x = x_values[min_index[0]]
    optimal_y = y_values[min_index[1]]

    return [optimal_x, optimal_y, 0], min_value
