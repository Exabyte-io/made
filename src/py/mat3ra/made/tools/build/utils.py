from typing import List

from mat3ra.made.material import Material


def merge_materials(materials: List[Material], distance_tolerance: float = 0.01) -> Material:
    """
    Merge materials into a single material.
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
    lattice = materials[0].lattice
    basis = []
    for material in materials:
        for atom in material.basis:
            basis.append(atom)

    name = "temp merged materials"
    new_config = {"name": name, "lattice": lattice, "basis": basis}
    return Material(new_config)
