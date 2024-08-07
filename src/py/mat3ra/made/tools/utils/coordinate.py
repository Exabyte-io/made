# Place all functions acting on coordinates
from typing import Callable, Dict, List, Tuple

import numpy as np


def is_coordinate_in_cylinder(
    coordinate: List[float], center_position: List[float], radius: float = 0.25, min_z: float = 0, max_z: float = 1
) -> bool:
    """
    Check if a coordinate is inside a cylinder.
    Args:
        coordinate (List[float]): The coordinate to check.
        center_position (List[float]): The coordinates of the center position.
        min_z (float): Lower limit of z-coordinate.
        max_z (float): Upper limit of z-coordinate.
        radius (float): The radius of the cylinder.

    Returns:
        bool: True if the coordinate is inside the cylinder, False otherwise.
    """
    return (coordinate[0] - center_position[0]) ** 2 + (coordinate[1] - center_position[1]) ** 2 <= radius**2 and (
        min_z <= coordinate[2] <= max_z
    )


def is_coordinate_in_sphere(coordinate: List[float], center_position: List[float], radius: float = 0.25) -> bool:
    """
    Check if a coordinate is inside a sphere.
    Args:
        coordinate (List[float]): The coordinate to check.
        center_position (List[float]): The coordinates of the center position.
        radius (float): The radius of the sphere.

    Returns:
        bool: True if the coordinate is inside the sphere, False otherwise.
    """
    np_coordinate = np.array(coordinate)
    np_center_position = np.array(center_position)
    distance_squared = np.sum((np_coordinate - np_center_position) ** 2)
    return distance_squared <= radius**2


def is_coordinate_in_box(
    coordinate: List[float], min_coordinate: List[float] = [0, 0, 0], max_coordinate: List[float] = [1, 1, 1]
) -> bool:
    """
    Check if a coordinate is inside a box.
    Args:
        coordinate (List[float]): The coordinate to check.
        min_coordinate (List[float]): The minimum coordinate of the box.
        max_coordinate (List[float]): The maximum coordinate of the box.
    Returns:
        bool: True if the coordinate is inside the box, False otherwise.
    """
    x_min, y_min, z_min = min_coordinate
    x_max, y_max, z_max = max_coordinate
    return x_min <= coordinate[0] <= x_max and y_min <= coordinate[1] <= y_max and z_min <= coordinate[2] <= z_max


def is_coordinate_within_layer(
    coordinate: List[float], center_position: List[float], direction_vector: List[float], layer_thickness: float
) -> bool:
    """
    Checks if a coordinate's projection along a specified direction vector
    is within a certain layer thickness centered around a given position.

    Args:
        coordinate (List[float]): The coordinate to check.
        center_position (List[float]): The coordinates of the center position.
        direction_vector (List[float]): The direction vector along which the layer thickness is defined.
        layer_thickness (float): The thickness of the layer along the direction vector.

    Returns:
        bool: True if the coordinate is within the layer thickness, False otherwise.
    """
    direction_norm = np.array(direction_vector) / np.linalg.norm(direction_vector)
    central_projection = np.dot(center_position, direction_norm)
    layer_thickness_frac = layer_thickness / np.linalg.norm(direction_vector)

    lower_bound = central_projection - layer_thickness_frac / 2
    upper_bound = central_projection + layer_thickness_frac / 2

    return lower_bound <= np.dot(coordinate, direction_norm) <= upper_bound


def is_coordinate_in_triangular_prism(
    coordinate: List[float],
    coordinate_1: List[float],
    coordinate_2: List[float],
    coordinate_3: List[float],
    min_z: float = 0,
    max_z: float = 1,
) -> bool:
    """
    Check if a coordinate is inside a triangular prism.
    Args:
        coordinate (List[float]): The coordinate to check.
        coordinate_1 (List[float]): The first coordinate of the triangle.
        coordinate_2 (List[float]): The second coordinate of the triangle.
        coordinate_3 (List[float]): The third coordinate of the triangle.
        min_z (float): Lower limit of z-coordinate.
        max_z (float): Upper limit of z-coordinate.

    Returns:
        bool: True if the coordinate is inside the triangular prism, False otherwise.
    """
    # convert to 3D coordinates at the origin XY plane
    coordinate_1.extend([0] * (3 - len(coordinate_1)))
    coordinate_2.extend([0] * (3 - len(coordinate_2)))
    coordinate_3.extend([0] * (3 - len(coordinate_3)))

    np_coordinate = np.array(coordinate)
    v1 = np.array(coordinate_1)
    v2 = np.array(coordinate_2)
    v3 = np.array(coordinate_3)

    v2_v1 = v2 - v1
    v3_v1 = v3 - v1
    coordinate_v1 = np_coordinate - v1

    # Compute dot products for the barycentric coordinates
    d00 = np.dot(v2_v1, v2_v1)
    d01 = np.dot(v2_v1, v3_v1)
    d11 = np.dot(v3_v1, v3_v1)
    d20 = np.dot(coordinate_v1, v2_v1)
    d21 = np.dot(coordinate_v1, v3_v1)

    # Calculate barycentric coordinates
    denom = d00 * d11 - d01 * d01
    v = (d11 * d20 - d01 * d21) / denom
    w = (d00 * d21 - d01 * d20) / denom
    u = 1.0 - v - w

    return (u >= 0) and (v >= 0) and (w >= 0) and (u + v + w <= 1) and (min_z <= np_coordinate[2] <= max_z)


def is_coordinate_behind_plane(
    coordinate: List[float], plane_normal: List[float], plane_point_coordinate: List[float]
) -> bool:
    """
    Check if a coordinate is behind a plane.
    Args:
        coordinate (List[float]): The coordinate to check.
        plane_normal (List[float]): The normal vector of the plane.
        plane_point_coordinate (List[float]): The coordinate of a point on the plane.

    Returns:
        bool: True if the coordinate is behind the plane, False otherwise.
    """
    np_coordinate = np.array(coordinate)
    np_plane_normal = np.array(plane_normal)
    np_plane_point = np.array(plane_point_coordinate)
    return np.dot(np_plane_normal, np_coordinate - np_plane_point) < 0


class CoordinateConditionBuilder:
    @staticmethod
    def create_condition(condition_type: str, evaluation_func: Callable, **kwargs) -> Tuple[Callable, Dict]:
        condition_json = {"type": condition_type, **kwargs}
        return lambda coordinate: evaluation_func(coordinate, **kwargs), condition_json

    @staticmethod
    def cylinder(center_position=None, radius: float = 0.25, min_z: float = 0, max_z: float = 1):
        if center_position is None:
            center_position = [0.5, 0.5]
        return CoordinateConditionBuilder.create_condition(
            condition_type="cylinder",
            evaluation_func=is_coordinate_in_cylinder,
            center_position=center_position,
            radius=radius,
            min_z=min_z,
            max_z=max_z,
        )

    @staticmethod
    def sphere(center_position=None, radius: float = 0.25):
        if center_position is None:
            center_position = [0.5, 0.5, 0.5]
        return CoordinateConditionBuilder.create_condition(
            condition_type="sphere",
            evaluation_func=is_coordinate_in_sphere,
            center_position=center_position,
            radius=radius,
        )

    @staticmethod
    def triangular_prism(
        position_on_surface_1: List[float] = [0, 0],
        position_on_surface_2: List[float] = [1, 0],
        position_on_surface_3: List[float] = [0, 1],
        min_z: float = 0,
        max_z: float = 1,
    ):
        return CoordinateConditionBuilder.create_condition(
            condition_type="prism",
            evaluation_func=is_coordinate_in_triangular_prism,
            coordinate_1=position_on_surface_1,
            coordinate_2=position_on_surface_2,
            coordinate_3=position_on_surface_3,
            min_z=min_z,
            max_z=max_z,
        )

    @staticmethod
    def box(min_coordinate=None, max_coordinate=None):
        if max_coordinate is None:
            max_coordinate = [1, 1, 1]
        if min_coordinate is None:
            min_coordinate = [0, 0, 0]
        return CoordinateConditionBuilder.create_condition(
            condition_type="box",
            evaluation_func=is_coordinate_in_box,
            min_coordinate=min_coordinate,
            max_coordinate=max_coordinate,
        )

    @staticmethod
    def plane(plane_normal: List[float], plane_point_coordinate: List[float]):
        return CoordinateConditionBuilder.create_condition(
            condition_type="plane",
            evaluation_func=is_coordinate_behind_plane,
            plane_normal=plane_normal,
            plane_point_coordinate=plane_point_coordinate,
        )
