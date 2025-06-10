from typing import List, Optional

import numpy as np
from mat3ra.code.array_with_ids import ArrayWithIds
from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from scipy.spatial import cKDTree

from mat3ra.made.basis import Basis, Coordinates
from mat3ra.made.material import Material
from mat3ra.made.tools.modify import translate_by_vector
from ...utils import AXIS_TO_INDEX_MAP


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
    if material1.lattice != material2.lattice and not merge_dangerously:
        raise ValueError("Lattices of the two materials must be the same.")
    merged_lattice = material1.lattice
    resolved_basis = merge_two_bases(material1.basis, material2.basis, distance_tolerance)

    name = material_name or "Merged Material"
    new_material = Material.create(
        {"name": name, "lattice": merged_lattice.to_dict(), "basis": resolved_basis.to_dict()}
    )
    # Handle None metadata gracefully
    metadata1 = material1.metadata or {}
    metadata2 = material2.metadata or {}
    new_material.metadata = {**metadata1, **metadata2}
    return new_material


def merge_materials(
    materials: List[Material],
    material_name: Optional[str] = None,
    distance_tolerance: float = 0.1,
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


def stack(materials: List[Material], direction: AxisEnum) -> Material:
    # use stack_two_materials
    return stack_two_materials(material_1=materials[0], material_2=materials[1], direction=direction)


def stack_two_materials(
    material_1: Material,
    material_2: Material,
    direction: AxisEnum,
) -> Material:
    """
    Stack two materials along a specified direction by expanding lattices and merging.

    Args:
        material_1: First material to stack.
        material_2: Second material to stack.
        direction: Direction along which to stack (x, y, or z axis).

    Returns:
        Material: Stacked material with combined lattice and atoms.
    """
    lattice_vector_index = AXIS_TO_INDEX_MAP[direction.value]

    material_1_lattice_vectors = material_1.lattice.vector_arrays
    material_2_lattice_vectors = material_2.lattice.vector_arrays

    stacked_lattice_vectors_values = [vec.copy() for vec in material_1_lattice_vectors]
    stacked_lattice_vectors_values[lattice_vector_index] = (
        np.array(material_1_lattice_vectors[lattice_vector_index])
        + np.array(material_2_lattice_vectors[lattice_vector_index])
    ).tolist()

    material_1_final_lattice_config = material_1.clone()
    material_1_final_lattice_config.set_new_lattice_vectors(
        lattice_vector1=stacked_lattice_vectors_values[0],
        lattice_vector2=stacked_lattice_vectors_values[1],
        lattice_vector3=stacked_lattice_vectors_values[2],
    )

    # Translate material2 so its atoms are positioned correctly relative to material1
    material2_translated = material_2.clone()
    # The translation amount is the original lattice vector of material1 in the stacking direction
    translation_vec = material_1_lattice_vectors[lattice_vector_index]
    material2_translated = translate_by_vector(material2_translated, translation_vec, use_cartesian_coordinates=True)

    stacked_material = merge_two_materials(
        material1=material_1_final_lattice_config,
        material2=material2_translated,
        distance_tolerance=0,
        merge_dangerously=True,  # Allows merging, assumes compatibility handled by stacking logic
        material_name=material_1.name,
    )

    return stacked_material
