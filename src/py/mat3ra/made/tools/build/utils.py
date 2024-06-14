from typing import List, Dict

import numpy as np
from mat3ra.made.material import Material
from mat3ra.made.utils import filter_array_with_id_value_by_ids

from mat3ra.made.tools.utils import convert_basis_to_crystal


def resolve_close_coordinates(basis: Dict, distance_tolerance: float = 0.01) -> Dict:
    """
    Resolve close coordinates in the basis by replacing them with the latest one.
    """
    elements = basis["elements"]
    coordinates = basis["coordinates"]
    labels = basis.get("labels", [])

    ids = set([coord["id"] for coord in coordinates])
    ids_to_remove = set()

    for i in range(1, len(coordinates)):
        for j in range(i):
            if (
                np.linalg.norm(np.array(coordinates[i]["value"]) - np.array(coordinates[j]["value"]))
                < distance_tolerance
            ):
                ids_to_remove.add(coordinates[j]["id"])

    ids_to_keep = list(ids - ids_to_remove)

    final_elements = filter_array_with_id_value_by_ids(elements, ids_to_keep)
    final_coordinates = filter_array_with_id_value_by_ids(coordinates, ids_to_keep)
    final_labels = filter_array_with_id_value_by_ids(labels, ids_to_keep)

    return {"elements": final_elements, "coordinates": final_coordinates, "labels": final_labels}


def merge_two_materials(material1: Material, material2: Material, distance_tolerance: float = 0.01) -> Material:
    """
    Merge two materials with the same lattice into a single material,
    replacing colliding atoms with the latest material's atoms.
    """

    material1 = material1.clone()
    material2 = material2.clone()
    if material1.lattice != material2.lattice:
        raise ValueError("Lattices of the two materials must be the same.")
    material1.basis = convert_basis_to_crystal(material1.basis)
    material2.basis = convert_basis_to_crystal(material2.basis)

    merged_lattice = material1.lattice
    merged_elements = material1.basis["elements"]
    merged_coordinates = material1.basis["coordinates"]
    merged_labels = material1.basis["labels"]

    for coordinate in material2.basis["coordinates"]:
        merged_coordinates.append(coordinate)

    for element in material2.basis["elements"]:
        merged_elements.append(element)

    if "labels" in material2.basis:
        for label in material2.basis["labels"]:
            merged_labels.append(label)

    merged_basis = {
        "elements": merged_elements,
        "coordinates": merged_coordinates,
        "labels": merged_labels,
    }

    resolved_basis = resolve_close_coordinates(merged_basis, distance_tolerance)

    merged_material_config = {
        "name": "Merged Material",
        "lattice": merged_lattice,
        "basis": resolved_basis,
    }
    return Material(merged_material_config)


def merge_materials(materials: List[Material], distance_tolerance: float = 0.01) -> Material:
    merged_material = materials[0]
    for material in materials[1:]:
        merged_material = merge_two_materials(merged_material, material, distance_tolerance)

    return merged_material
