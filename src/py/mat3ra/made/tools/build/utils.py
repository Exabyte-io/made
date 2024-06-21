from typing import List

import numpy as np
from mat3ra.made.basis.basis import Basis
from mat3ra.made.material import Material

from mat3ra.made.tools.utils import convert_basis_to_crystal


def resolve_close_coordinates_basis(basis: Basis, distance_tolerance: float = 0.01) -> Basis:
    """
    Find all atoms that are within distance tolerance and only keep the last one, remove other sites
    """
    coordinates = basis.coordinates.to_array_of_values_with_ids()
    ids = set([coord.id for coord in coordinates])
    ids_to_remove = set()

    for i in range(1, len(coordinates)):
        for j in range(i):
            if np.linalg.norm(np.array(coordinates[i].value) - np.array(coordinates[j].value)) < distance_tolerance:
                ids_to_remove.add(coordinates[j].id)

    ids_to_keep = list(ids - ids_to_remove)
    basis.filter_atoms_by_ids(ids_to_keep)
    return basis


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
    merged_elements = material1.basis.elements
    merged_coordinates = material1.basis.coordinates
    merged_labels = material1.basis.labels

    for coordinate in material2.basis.coordinates.values:
        merged_coordinates.add_item(coordinate)

    for element in material2.basis.elements.values:
        merged_elements.add_item(element)

    if material2.basis.labels:
        for label in material2.basis.labels.values:
            merged_labels.add_item(label)

    merged_basis = Basis(
        elements=merged_elements,
        coordinates=merged_coordinates,
        units=material1.basis.units,
        cell=material1.basis.cell,
        labels=merged_labels,
    )
    resolved_basis = resolve_close_coordinates_basis(merged_basis, distance_tolerance)

    name = "Merged Material"
    new_material = Material.create(
        {"name": name, "lattice": merged_lattice.to_json(), "basis": resolved_basis.to_json()}
    )
    return new_material


def merge_materials(materials: List[Material], distance_tolerance: float = 0.01) -> Material:
    merged_material = materials[0]
    for material in materials[1:]:
        merged_material = merge_two_materials(merged_material, material, distance_tolerance)

    return merged_material
