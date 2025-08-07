from typing import List, Tuple, Union

import numpy as np

from mat3ra.made.material import Material
from .enums import StrainModes
from .....build_components import MaterialWithBuildMetadata
from .....convert.interface_parts_enum import INTERFACE_LABELS_MAP
from .....modify import interface_displace_part, filter_by_label
from .....optimize import evaluate_calculator_on_xy_grid
from .....third_party import PymatgenInterface


def remove_duplicate_interfaces(
    interfaces: List[PymatgenInterface], strain_mode: StrainModes = StrainModes.mean_abs_strain
):
    def are_interfaces_duplicate(interface1: PymatgenInterface, interface2: PymatgenInterface):
        are_sizes_equivalent = interface1.num_sites == interface2.num_sites and np.allclose(
            interface1.interface_properties[strain_mode], interface2.interface_properties[strain_mode]
        )
        are_strains_equivalent = np.allclose(
            interface1.interface_properties[strain_mode], interface2.interface_properties[strain_mode]
        )
        return are_sizes_equivalent and are_strains_equivalent

    filtered_interfaces = [interfaces[0]] if interfaces else []

    for interface in interfaces[1:]:
        if not any(are_interfaces_duplicate(interface, unique_interface) for unique_interface in filtered_interfaces):
            filtered_interfaces.append(interface)
    return filtered_interfaces


def get_slab(interface: Material, part: str = "film"):
    try:
        return filter_by_label(interface, INTERFACE_LABELS_MAP[part])
    except ValueError:
        raise ValueError(f"Material does not contain label for {part}.")


def get_optimal_film_displacement(
    material: Union[Material, MaterialWithBuildMetadata],
    calculator,
    grid_size_xy: Tuple[int, int] = (10, 10),
    grid_offset_position: List[float] = [0, 0],
    grid_range_x=(-0.5, 0.5),
    grid_range_y=(-0.5, 0.5),
    use_cartesian_coordinates=False,
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
        calculator: The calculator to use for the calculation of the interaction energy.

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
