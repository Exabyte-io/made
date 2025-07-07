from typing import List

import numpy as np
from pydantic import BaseModel


class FunctionHolder(BaseModel):
    def apply_function(self, coordinate: List[float]) -> float:
        """
        Get the value of the function at the given coordinate.
        """
        raise NotImplementedError

    def calculate_derivative(self, coordinate: List[float], axis: str) -> float:
        """
        Get the derivative of the function at the given coordinate
        """
        raise NotImplementedError

    def calculate_arc_length(self, a: float, b: float) -> float:
        """
        Get the arc length of the function between a and b.
        """
        raise NotImplementedError

    def get_json(self) -> dict:
        """
        Get the json representation of the function holder.
        """
        raise NotImplementedError


def calculate_required_vacuum_for_rotation(material, angle_deg, vacuum_x, vacuum_y):
    """
    Calculate the minimum cell size in x and y so that, after rotation by angle_deg,
    the farthest atom will have at least vacuum_x/y from the cell edge.
    Returns (required_length_x, required_length_y, center_shift)
    """
    coords = np.array(material.coordinates_array)
    # Convert fractional to cartesian
    lattice = np.array(material.lattice.vector_arrays)
    cart_coords = np.dot(coords, lattice)
    # Rotation matrix (about z)
    theta = np.radians(angle_deg)
    rot_matrix = np.array(
        [
            [np.cos(theta), -np.sin(theta), 0],
            [np.sin(theta), np.cos(theta), 0],
            [0, 0, 1],
        ]
    )
    rotated = np.dot(cart_coords, rot_matrix.T)
    min_x, max_x = np.min(rotated[:, 0]), np.max(rotated[:, 0])
    min_y, max_y = np.min(rotated[:, 1]), np.max(rotated[:, 1])
    length_x = (max_x - min_x) + 2 * vacuum_x
    length_y = (max_y - min_y) + 2 * vacuum_y
    # Center shift to move structure to center of new cell
    center_x = (max_x + min_x) / 2
    center_y = (max_y + min_y) / 2
    center_shift = np.array([length_x / 2 - center_x, length_y / 2 - center_y, 0])
    return length_x, length_y, center_shift
