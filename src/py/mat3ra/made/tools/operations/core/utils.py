from typing import List, Optional

import numpy as np
from mat3ra.code.array_with_ids import ArrayWithIds
from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from mat3ra.made.basis import Basis, Coordinates
from mat3ra.made.material import Material

from ...convert import from_ase, to_ase
from ...third_party import ase_make_supercell
from ...utils import decorator_convert_supercell_matrix_2x2_to_3x3


def merge_two_bases(basis1: Basis, basis2: Basis, distance_tolerance: float) -> Basis:
    basis1.to_crystal()
    basis2.to_crystal()

    merged_elements_values = basis1.elements.values + basis2.elements.values
    merged_coordinates_values = basis1.coordinates.values + basis2.coordinates.values
    merged_labels_values = basis1.labels.values + basis2.labels.values if basis1.labels and basis2.labels else []
    merged_constraints_values = basis1.constraints.values + basis2.constraints.values

    new_basis = basis1.clone()
    new_basis.elements = ArrayWithIds.from_values(values=merged_elements_values)
    new_basis.coordinates = Coordinates.from_values(values=merged_coordinates_values)
    new_basis.labels = ArrayWithIds.from_values(values=merged_labels_values)
    new_basis.constraints = ArrayWithIds.from_values(values=merged_constraints_values)
    new_basis.resolve_colliding_coordinates(tolerance=distance_tolerance)

    return new_basis


def merge_two_materials(
    material1: Material,
    material2: Material,
    material_name: Optional[str] = None,
    distance_tolerance: float = 0.1,
    merge_dangerously=False,
) -> Material:
    if not np.allclose(material1.lattice.vector_arrays, material2.lattice.vector_arrays) and not merge_dangerously:
        raise ValueError("Lattices of the two materials must be the same.")
    merged_lattice = material1.lattice
    resolved_basis = merge_two_bases(material1.basis, material2.basis, distance_tolerance)

    name = material_name or "Merged Material"
    new_material = Material.create(
        {"name": name, "lattice": merged_lattice.to_dict(), "basis": resolved_basis.to_dict()}
    )
    return new_material


def should_skip_stacking(
    material_1: Material,
    material_2: Material,
    lattice_vector_index: int,
) -> bool:
    """
    Determine if stacking should be skipped due to zero thickness in material_2.

    Also validates that materials are compatible for stacking.

    Args:
        material_1: First material to stack.
        material_2: Second material to stack.
        lattice_vector_index: Index of the lattice vector in the stacking direction.

    Returns:
        bool: True if stacking should be skipped (material_2 has zero thickness).

    Raises:
        ValueError: If in-plane lattice vectors are not the same.
    """
    material_1_lattice_vectors = material_1.lattice.vector_arrays
    material_2_lattice_vectors = material_2.lattice.vector_arrays

    is_zero_thickness = np.allclose(material_2_lattice_vectors[lattice_vector_index], [0, 0, 0])
    if is_zero_thickness:
        return True

    # Check that in-plane vectors are the same
    all_indices = [0, 1, 2]
    in_plane_indices = [i for i in all_indices if i != lattice_vector_index]

    material_1_in_plane_vectors = np.array(material_1_lattice_vectors)[in_plane_indices]
    material_2_in_plane_vectors = np.array(material_2_lattice_vectors)[in_plane_indices]

    if not np.allclose(material_1_in_plane_vectors, material_2_in_plane_vectors):
        raise ValueError("In-plane lattice vectors of the two materials must be the same for stacking.")

    return False


@decorator_convert_supercell_matrix_2x2_to_3x3
def mirror(material: Material, direction: AxisEnum = AxisEnum.z) -> Material:
    """
    Mirrors the material along the specified axis by applying a right-handed supercell transformation.
    """
    supercell_matrix = {
        AxisEnum.x: [[-1, 0, 0], [0, 0, 1], [0, 1, 0]],
        AxisEnum.y: [[0, 0, 1], [0, -1, 0], [1, 0, 0]],
        AxisEnum.z: [[0, 1, 0], [1, 0, 0], [0, 0, -1]],
    }.get(direction)

    atoms = to_ase(material)

    supercell_atoms = ase_make_supercell(atoms, supercell_matrix)
    new_material = Material.create(from_ase(supercell_atoms))
    if material.metadata:
        new_material.metadata = material.metadata
    new_material.name = material.name
    return new_material
