from typing import List, Optional, Union

import numpy as np
from mat3ra.code.array_with_ids import ArrayWithIds
from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from scipy.spatial import cKDTree

from mat3ra.made.basis import Basis, Coordinates
from mat3ra.made.material import Material
from mat3ra.made.tools.modify import translate_by_vector
from .slab.configuration import VacuumConfiguration
from .supercell import create_supercell
from ..modify import filter_by_box
from ..utils import AXIS_TO_INDEX_MAP


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
    merged_constraints_values = basis1.constraints.values + basis2.constraints.values

    new_basis = basis1.clone()
    new_basis.elements = ArrayWithIds.from_values(values=merged_elements_values)
    new_basis.coordinates = Coordinates.from_values(values=merged_coordinates_values)
    new_basis.labels = ArrayWithIds.from_values(values=merged_labels_values)
    new_basis.constraints = ArrayWithIds.from_values(values=merged_constraints_values)

    resolved_basis = resolve_close_coordinates_basis(new_basis, distance_tolerance)

    return resolved_basis


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


def double_and_filter_material(material: Material, start: List[float], end: List[float]) -> Material:
    """
    Double the material and filter it by a box defined by the start and end coordinates.
    Args:
        material (Material): The material to double and filter.
        start (List[float]): The start coordinates of the box.
        end (List[float]): The end coordinates of the box.
    Returns:
        Material: The filtered material.
    """
    material_doubled = create_supercell(material, scaling_factor=[2, 1, 1])
    return filter_by_box(material_doubled, start, end)


def expand_lattice_vectors(material: Material, gap: float, direction: int = 0) -> Material:
    """
    Expand the lattice vectors of the material in the specified direction by the given gap.

    Args:
        material (Material): The material whose lattice vectors are to be expanded.
        gap (float): The gap by which to expand the lattice vector.
        direction (int): The index of the lattice vector to expand (0, 1, or 2).
    """
    new_lattice_vectors = material.lattice.vector_arrays
    new_lattice_vectors[direction][direction] += gap
    material.set_new_lattice_vectors(
        lattice_vector1=new_lattice_vectors[0],
        lattice_vector2=new_lattice_vectors[1],
        lattice_vector3=new_lattice_vectors[2],
    )
    return material


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
    )

    return stacked_material


def stack_two_materials_xy(
    phase_1_material: Material,
    phase_2_material: Material,
    gap: float,
    edge_inclusion_tolerance: Optional[float] = 1.0,
    distance_tolerance: float = 1.0,
) -> Material:
    """
    Stack two materials laterally with translation along x-axis with a gap between them.

    Works correctly only for materials with the same lattice vectors (commensurate lattices).
     Args:
        phase_1_material (Material): The first material.
        phase_2_material (Material): The second material.
        gap (float): The gap between the two materials, in angstroms.
        edge_inclusion_tolerance (float): The tolerance to include atoms on the edge of the phase, in angstroms.
        distance_tolerance (float): The distance tolerance to remove atoms that are too close, in angstroms.

    Returns:
        Material: The merged material.
    """
    edge_inclusion_tolerance_crystal = abs(
        phase_1_material.basis.cell.convert_point_to_crystal([edge_inclusion_tolerance, 0, 0])[0]
    )

    phase_1_material = double_and_filter_material(
        phase_1_material, [0 - edge_inclusion_tolerance_crystal, 0, 0], [0.5 + edge_inclusion_tolerance_crystal, 1, 1]
    )

    phase_2_material = double_and_filter_material(
        phase_2_material, [0.5 - edge_inclusion_tolerance_crystal, 0, 0], [1 + edge_inclusion_tolerance_crystal, 1, 1]
    )

    phase_1_material = expand_lattice_vectors(phase_1_material, gap)
    phase_2_material = expand_lattice_vectors(phase_2_material, gap)

    phase_2_material = translate_by_vector(phase_2_material, [gap / 2, 0, 0], use_cartesian_coordinates=True)
    interface = merge_materials(
        [phase_1_material, phase_2_material], distance_tolerance=distance_tolerance, merge_dangerously=True
    )
    return interface


def stack_two_components(
    component1: Union[Material, "VacuumConfiguration"],
    component2: Union[Material, "VacuumConfiguration"],
    direction: AxisEnum,
) -> Material:
    from .slab.configuration import VacuumConfiguration

    if isinstance(component1, Material):
        reference_material = component1
    elif isinstance(component2, Material):
        reference_material = component2
    else:
        raise ValueError("At least one component must be a Material to serve as reference for configurations")

    # Convert components to materials if they are configurations
    if isinstance(component1, VacuumConfiguration):
        material1 = create_vacuum_material(reference_material, component1)
    else:
        material1 = component1

    if isinstance(component2, VacuumConfiguration):
        material2 = create_vacuum_material(reference_material, component2)
    else:
        material2 = component2

    return stack_two_materials(material1, material2, direction)
