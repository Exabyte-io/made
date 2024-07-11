from typing import List, Optional

from mat3ra.made.basis import Basis
from mat3ra.made.material import Material

from ..utils import get_distance_between_coordinates


def resolve_close_coordinates_basis(basis: Basis, distance_tolerance: float = 0.01) -> Basis:
    """
    Find all atoms that are within distance tolerance and only keep the last one, remove other sites
    """
    coordinates = basis.coordinates.to_array_of_values_with_ids()
    ids = set(basis.coordinates.ids)
    ids_to_remove = set()

    for i in range(1, len(coordinates)):
        for j in range(i):
            if get_distance_between_coordinates(coordinates[i].value, coordinates[j].value) < distance_tolerance:
                ids_to_remove.add(coordinates[j].id)

    ids_to_keep = list(ids - ids_to_remove)
    basis.filter_atoms_by_ids(ids_to_keep)
    return basis


def merge_two_bases(basis1: Basis, basis2: Basis, distance_tolerance: float) -> Basis:
    basis1.to_crystal()
    basis2.to_crystal()

    merged_elements = basis1.elements
    merged_coordinates = basis1.coordinates
    merged_labels = basis1.labels

    for coordinate in basis2.coordinates.values:
        merged_coordinates.add_item(coordinate)

    for element in basis2.elements.values:
        merged_elements.add_item(element)

    if basis2.labels:
        for label in basis2.labels.values:
            merged_labels.add_item(label)

    merged_basis = Basis(
        elements=merged_elements,
        coordinates=merged_coordinates,
        units=basis1.units,
        cell=basis1.cell,
        labels=merged_labels,
    )
    resolved_basis = resolve_close_coordinates_basis(merged_basis, distance_tolerance)

    return resolved_basis


def merge_two_materials(
    material1: Material,
    material2: Material,
    material_name: Optional[str],
    distance_tolerance: float,
    merge_dangerously=False,
) -> Material:
    """
    Merge two materials with the same lattice into a single material,
    replacing colliding atoms with the latest material's atoms.

    Args:
        material1 (Material): First material to merge.
        material2 (Material): Second material to merge.
        material_name (str): Name of the merged material.
        distance_tolerance (float): Distance tolerance for merging atoms.
        merge_dangerously (bool): If True, merge materials with different lattices.
    """

    material1 = material1.clone()
    material2 = material2.clone()
    if material1.lattice != material2.lattice and not merge_dangerously:
        raise ValueError("Lattices of the two materials must be the same.")
    merged_lattice = material1.lattice
    resolved_basis = merge_two_bases(material1.basis, material2.basis, distance_tolerance)

    name = material_name or "Merged Material"
    new_material = Material.create(
        {"name": name, "lattice": merged_lattice.to_json(), "basis": resolved_basis.to_json()}
    )
    new_material.metadata = {**material1.metadata, **material2.metadata}
    return new_material


def merge_materials(
    materials: List[Material],
    material_name: Optional[str] = None,
    distance_tolerance: float = 0.01,
    merge_dangerously=False,
) -> Material:
    merged_material = materials[0]
    for material in materials[1:]:
        merged_material = merge_two_materials(
            merged_material, material, material_name, distance_tolerance, merge_dangerously
        )

    return merged_material
