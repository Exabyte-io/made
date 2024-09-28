import numpy as np
from scipy.spatial import cKDTree
from typing import List, Optional
from mat3ra.made.basis import Basis
from mat3ra.made.material import Material

from .supercell import create_supercell
from ..modify import filter_by_box, translate_by_vector
from ...utils import ArrayWithIds


def resolve_close_coordinates_basis(basis: Basis, distance_tolerance: float = 0.1) -> Basis:
    """
    Find all atoms that are within distance tolerance and only keep the last one, remove other sites

    Args:
        basis (Basis): The basis to resolve.
        distance_tolerance (float): The distance tolerance in angstroms.
    """
    basis.to_cartesian()
    coordinates = np.array(basis.coordinates.values)
    tree = cKDTree(coordinates)
    colliding_pairs = tree.query_pairs(r=distance_tolerance)

    ids_to_remove = set()
    if len(colliding_pairs) == 0:
        basis.to_crystal()
        return basis
    for index_1, index_2 in colliding_pairs:
        ids_to_remove.add(index_1)

    ids_to_keep = [id for id in basis.coordinates.ids if id not in ids_to_remove]
    basis = basis.filter_atoms_by_ids(ids_to_keep)
    basis.to_crystal()
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
    """
    Merge multiple materials into a single material.

    If some of the atoms are considered too close within a tolerance, only the last atom is kept.

    Args:
        materials (List[Material]): List of materials to merge.
        material_name (Optional[str]): Name of the merged material.
        distance_tolerance (float): The tolerance to replace atoms that are considered too close with respect
            to the coordinates in the last material in the list, in angstroms.
        merge_dangerously (bool): If True, the lattices are merged "as is" with no sanity checks.
    Returns:
        Material: The merged material.
    """
    merged_material = materials[0]
    for material in materials[1:]:
        merged_material = merge_two_materials(
            merged_material, material, material_name, distance_tolerance, merge_dangerously
        )

    return merged_material


def merge_two_materials_laterally(
    phase_1_material: Material,
    phase_2_material: Material,
    gap: float,
    cut_tolerance: Optional[float] = None,
    distance_tolerance: float = 1.0,
) -> Material:
    """
    Merge two materials laterally with translation along x axis with a gap between them.
    Args:
        phase_1_material (Material): The first material.
        phase_2_material (Material): The second material.
        gap (float): The gap between the two materials, in angstroms.
        cut_tolerance (float): The cut tolerance to include atoms on the edge of cut, in crystal coordinates.
        distance_tolerance (float): The distance tolerance to remove atoms that are too close, in angstroms.

    Returns:
        Material: The merged material.
    """
    # Default cut tolerance is half of the gap to allow atoms to fill the gap if they are on the edge
    cut_tolerance = cut_tolerance or phase_1_material.basis.cell.convert_point_to_cartesian([gap, 0, 0])[2] / 2

    phase_1_material_doubled = create_supercell(phase_1_material, scaling_factor=[2, 1, 1])
    phase_1_material = filter_by_box(phase_1_material_doubled, [0 - cut_tolerance, 0, 0], [0.5, 1, 1])

    phase_2_material_doubled = create_supercell(phase_2_material, scaling_factor=[2, 1, 1])
    phase_2_material = filter_by_box(phase_2_material_doubled, [0.5 - cut_tolerance, 0, 0], [1, 1, 1])

    new_lattice_vectors_1 = phase_1_material.lattice.vector_arrays
    new_lattice_vectors_1[0][0] += gap

    new_lattice_vectors_2 = phase_2_material.lattice.vector_arrays
    new_lattice_vectors_2[0][0] += gap

    phase_1_material.set_new_lattice_vectors(
        lattice_vector1=new_lattice_vectors_1[0],
        lattice_vector2=new_lattice_vectors_1[1],
        lattice_vector3=new_lattice_vectors_1[2],
    )
    phase_2_material.set_new_lattice_vectors(
        lattice_vector1=new_lattice_vectors_2[0],
        lattice_vector2=new_lattice_vectors_2[1],
        lattice_vector3=new_lattice_vectors_2[2],
    )

    phase_2_material = translate_by_vector(phase_2_material, [gap / 2, 0, 0], use_cartesian_coordinates=True)
    interface = merge_materials([phase_1_material, phase_2_material], distance_tolerance=distance_tolerance)
    return interface
