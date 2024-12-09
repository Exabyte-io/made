from functools import wraps
from typing import Union

import numpy as np
from mat3ra.made.material import Material
from mat3ra.made.utils import ArrayWithIds

from .material import MaterialWithCrystalSites


def decorator_handle_periodic_boundary_conditions(cutoff):
    """
    Decorator to handle periodic boundary conditions.

    Copies atoms near boundaries within the cutoff distance beyond the opposite side of the cell
    creating the effect of periodic boundary conditions for edge atoms.

    Results of the function are filtered to remove atoms or coordinates outside the original cell.

    Args:
        cutoff (float): The cutoff distance for a border slice in crystal coordinates.

    Returns:
        Callable: The decorated function.
    """

    def decorator(func):
        @wraps(func)
        def wrapper(material, *args, **kwargs):
            augmented_material, last_id = augment_material_with_periodic_images(material, cutoff)
            result = func(augmented_material, *args, **kwargs)
            original_ids = material.basis.coordinates.ids
            if isinstance(result, ArrayWithIds):
                result.filter_by_ids(original_ids)
            if isinstance(result, list):
                if all(isinstance(x, int) for x in result):
                    result = [id for id in result if id <= last_id]
                elif all(isinstance(x, list) and len(x) == 3 for x in result):
                    result = [coord for coord in result if all(0 <= c < 1 for c in coord)]
            return result

        return wrapper

    return decorator


def augment_material_with_periodic_images(material: Union[Material, MaterialWithCrystalSites], cutoff: float = 1.0):
    """
    Augment the material's dataset by adding atoms from periodic images within a cutoff distance from the boundaries by
    copying them to the opposite side of the cell, translated by the cell vector beyond the boundary.

    Args:
        material (Material): The material to augment.
        cutoff (float): The cutoff value for filtering atoms near boundaries, in crystal coordinates.

    Returns:
        Tuple[Material, int]: The augmented material and the original count of atoms.
    """
    last_id = material.basis.coordinates.ids[-1]
    coordinates = np.array(material.basis.coordinates.values)
    elements = np.array(material.basis.elements.values)
    augmented_material = material.clone()
    new_basis = augmented_material.basis.copy()

    for axis in range(3):
        for direction in [-1, 1]:
            mask = (coordinates[:, axis] < cutoff) if direction == 1 else (coordinates[:, axis] > (1 - cutoff))
            filtered_coordinates = coordinates[mask]
            filtered_elements = elements[mask]
            translation_vector = np.zeros(3)
            translation_vector[axis] = direction
            translated_coordinates = filtered_coordinates + translation_vector
            for coord, elem in zip(translated_coordinates, filtered_elements):
                new_basis.add_atom(elem, coord)

    augmented_material.basis = new_basis
    return augmented_material, last_id