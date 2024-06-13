from typing import List, Dict
from itertools import zip_longest

import numpy as np
from mat3ra.made.material import Material


def merge_materials(materials: List[Material], distance_tolerance: float = 0.01) -> Material:
    """
    Merge materials into a single material using first material's lattice.
    Any collisions in coordinates are resolved by keeping the last atom.
    Args:
        materials: List of materials to merge.
        distance_tolerance: Tolerance for merging atoms in Angstrom.

    Returns:
        Material: Merged material.
    """
    # lattice from first
    # basis from all
    # remove duplicated atoms by coordinate within tolerance
    # prioritize last atom in case of duplicates

    # Get lattice from first material
    merged_lattice = materials[0].lattice
    merged_basis: Dict = {"elements": [], "coordinates": [], "labels": []}
    existing_coords: List = []

    for material in materials:
        # Use zip_longest to handle different lengths, filling missing labels with None
        for element, coord, label in zip_longest(
            material.basis["elements"], material.basis["coordinates"], material.basis["labels"], fillvalue=None
        ):
            new_pos = np.array(coord["value"])  # type: ignore
            is_duplicate = False

            # Check if this new position clashes with existing positions
            for i, existing in enumerate(existing_coords):
                if np.linalg.norm(new_pos - np.array(existing["value"])) < distance_tolerance:
                    # If a duplicate is found, replace the existing atom
                    merged_basis["elements"][i] = element
                    merged_basis["coordinates"][i] = coord
                    merged_basis["labels"][i] = label if label is not None else merged_basis["labels"][i]
                    is_duplicate = True
                    break

            if not is_duplicate:
                # Add new atom if no duplicate was found
                merged_basis["elements"].append(element)
                merged_basis["coordinates"].append(coord)
                merged_basis["labels"].append(label if label is not None else "")
                existing_coords.append(coord)

    merged_material_config = {"name": "Merged Material", "lattice": merged_lattice, "basis": merged_basis}
    return Material(merged_material_config)
