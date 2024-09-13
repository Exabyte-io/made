from typing import Union, List, Optional, Tuple, Callable

import numpy as np

from mat3ra.made.material import Material
from ...calculate import calculate_film_substrate_interaction_metric
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
from ...utils import get_sum_of_inverse_distances_squared


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
    calculator: Callable = calculate_film_substrate_interaction_metric,
):
    x_values, y_values, results_matrix = calculate_on_xy_grid(
        material,
        modifier=displace_interface_part,
        modifier_parameters={},
        calculator=calculator,
        calculator_parameters={"shadowing_radius": 2.5, "metric_function": get_sum_of_inverse_distances_squared},
        grid_size_xy=grid_size_xy,
        grid_offset_position=grid_offset_position,
        grid_range_x=grid_range_x,
        grid_range_y=grid_range_y,
        use_cartesian_coordinates=use_cartesian_coordinates,
    )
    min_index = np.unravel_index(np.argmin(results_matrix), results_matrix.shape)

    optimal_x = x_values[min_index[0]]
    optimal_y = y_values[min_index[1]]

    return [optimal_x, optimal_y, 0]
