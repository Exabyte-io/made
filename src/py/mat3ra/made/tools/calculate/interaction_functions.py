import numpy as np


def sum_of_inverse_distances_squared(
    coordinates_1: np.ndarray, coordinates_2: np.ndarray, epsilon: float = 1e-12
) -> float:
    """
    Calculate the sum of inverse squares of distances between two sets of coordinates.

    Args:
        coordinates_1 (np.ndarray): The first set of coordinates, shape (N1, 3).
        coordinates_2 (np.ndarray): The second set of coordinates, shape (N2, 3).
        epsilon (float): Small value to prevent division by zero.

    Returns:
        float: The calculated sum.
    """
    differences = coordinates_1[:, np.newaxis, :] - coordinates_2[np.newaxis, :, :]  # Shape: (N1, N2, 3)
    distances_squared = np.sum(differences**2, axis=2)  # Shape: (N1, N2)
    distances_squared = np.where(distances_squared == 0, epsilon, distances_squared)
    inv_distances_squared = -1 / distances_squared
    total = np.sum(inv_distances_squared)
    return float(total)
