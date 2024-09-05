import numpy as np
from scipy.spatial import cKDTree
from typing import List, Optional
from mat3ra.made.basis import Basis
from mat3ra.made.material import Material
from ...utils import ArrayWithIds


def resolve_close_coordinates_basis(basis: Basis, distance_tolerance: float = 0.01) -> Basis:
    """
    Find all atoms that are within distance tolerance and only keep the last one, remove other sites
    """
    coordinates = np.array(basis.coordinates.values)
    tree = cKDTree(coordinates)
    colliding_pairs = tree.query_pairs(r=distance_tolerance)

    ids_to_remove = set()
    if len(colliding_pairs) == 0:
        return basis
    for index_1, index_2 in colliding_pairs:
        ids_to_remove.add(index_1)

    ids_to_keep = [id for id in basis.coordinates.ids if id not in ids_to_remove]
    basis = basis.filter_atoms_by_ids(ids_to_keep)
    return basis


def merge_two_bases(basis1: Basis, basis2: Basis, distance_tolerance: float) -> Basis:
    basis1.to_crystal()
    basis2.to_crystal()

    merged_elements_values = basis1.elements.values + basis2.elements.values
    merged_coordinates_values = basis1.coordinates.values + basis2.coordinates.values
    merged_labels_values = basis1.labels.values + basis2.labels.values if basis1.labels and basis2.labels else []

    merged_basis = Basis(
        elements=ArrayWithIds.from_values(values=merged_elements_values),
        coordinates=ArrayWithIds.from_values(values=merged_coordinates_values),
        units=basis1.units,
        cell=basis1.cell,
        labels=ArrayWithIds.from_values(values=merged_labels_values),
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
