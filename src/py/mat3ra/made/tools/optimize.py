from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np
from mat3ra.made.material import Material


def evaluate_calculator_on_xy_grid(
    material: Material,
    calculator_function: Callable[[Material], Any],
    modifier: Optional[Callable] = None,
    modifier_parameters: Dict[str, Any] = {},
    grid_size_xy: Tuple[int, int] = (10, 10),
    grid_offset_position: List[float] = [0, 0],
    grid_range_x: Tuple[float, float] = (-0.5, 0.5),
    grid_range_y: Tuple[float, float] = (-0.5, 0.5),
    use_cartesian_coordinates: bool = False,
) -> Tuple[List[np.ndarray], np.ndarray]:
    """
    Calculate a property on a grid of x-y positions.

    Args:
        material (Material): The material object.
        modifier (Callable): The modifier function to apply to the material.
        modifier_parameters (Dict[str, Any]): The parameters to pass to the modifier.
        calculator_function (Callable): The calculator function to apply to the modified material.
        grid_size_xy (Tuple[int, int]): The size of the grid in x and y directions.
        grid_offset_position (List[float]): The offset position of the grid, in Angstroms or crystal coordinates.
        grid_range_x (Tuple[float, float]): The range to search in x direction, in Angstroms or crystal coordinates.
        grid_range_y (Tuple[float, float]): The range to search in y direction, in Angstroms or crystal coordinates.
        use_cartesian_coordinates (bool): Whether to use Cartesian coordinates.

    Returns:
        Tuple[List[np.ndarray[float]], np.ndarray[float]]: The x-y positions and the calculated property values.
    """
    x_values = np.linspace(grid_range_x[0], grid_range_x[1], grid_size_xy[0]) + grid_offset_position[0]
    y_values = np.linspace(grid_range_y[0], grid_range_y[1], grid_size_xy[1]) + grid_offset_position[1]

    xy_matrix = [x_values, y_values]
    results_matrix = np.zeros(grid_size_xy)

    for i, x in enumerate(x_values):
        for j, y in enumerate(y_values):
            if modifier is None:
                modified_material = material
            else:
                modified_material = modifier(
                    material,
                    displacement=[x, y, 0],
                    use_cartesian_coordinates=use_cartesian_coordinates,
                    **modifier_parameters,
                )
            result = calculator_function(modified_material)
            results_matrix[i, j] = result

    return xy_matrix, results_matrix
