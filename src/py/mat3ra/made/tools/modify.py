from enum import Enum
from typing import Callable, List, Literal, Optional, Union

import numpy as np
from mat3ra.made.material import Material
from mat3ra.made.tools.build.supercell import create_supercell
from mat3ra.made.tools.utils import transform_coordinate_to_supercell

from .analyze import (
    get_atom_indices_with_condition_on_coordinates,
    get_atom_indices_within_radius_pbc,
    get_atomic_coordinates_extremum,
    get_undercoordinated_atoms,
)
from .convert import from_ase, to_ase
from .third_party import ase_add_vacuum
from .utils.coordinate import (
    is_coordinate_in_box,
    is_coordinate_in_cylinder,
    is_coordinate_in_triangular_prism,
    is_coordinate_within_layer,
)


def filter_by_label(material: Material, label: Union[int, str]) -> Material:
    """
    Filter out only atoms corresponding to the label.

    Args:
        material (Material): The material object to filter.
        label (int|str): The tag/label to filter by.

    Returns:
        Material: The filtered material object.
    """
    new_material = material.clone()
    labels_array = new_material.basis.labels.to_array_of_values_with_ids()
    filtered_label_ids = [_label.id for _label in labels_array if _label.value == label]
    new_basis = new_material.basis
    new_basis.filter_atoms_by_ids(filtered_label_ids)
    new_material.basis = new_basis
    return new_material


def translate_to_z_level(
    material: Material, z_level: Optional[Literal["top", "bottom", "center"]] = "bottom"
) -> Material:
    """
    Translate atoms to the specified z-level.

    Args:
        material (Material): The material object to normalize.
        z_level (str): The z-level to translate the atoms to (top, bottom, center)
    Returns:
        Material: The translated material object.
    """
    min_z = get_atomic_coordinates_extremum(material, "min")
    max_z = get_atomic_coordinates_extremum(material)
    if z_level == "top":
        material = translate_by_vector(material, vector=[0, 0, 1 - max_z])
    elif z_level == "bottom":
        material = translate_by_vector(material, vector=[0, 0, -min_z])
    elif z_level == "center":
        material = translate_by_vector(material, vector=[0, 0, (1 - min_z - max_z) / 2])
    return material


def translate_by_vector(
    material: Material,
    vector: Optional[List[float]] = None,
    use_cartesian_coordinates: bool = False,
) -> Material:
    """
    Translate atoms by a vector.

    Args:
        material (Material): The material object to normalize.
        vector (List[float]): The vector to translate the atoms by (in crystal coordinates by default).
        use_cartesian_coordinates (bool): Whether to use cartesian coordinates.
    Returns:
        Material: The translated material object.
    """
    if not use_cartesian_coordinates:
        vector = material.basis.cell.convert_point_to_cartesian(vector)

    if vector is None:
        vector = [0, 0, 0]

    atoms = to_ase(material)
    # ASE accepts cartesian coordinates for translation
    atoms.translate(tuple(vector))
    return Material(from_ase(atoms))


def translate_to_center(material: Material, axes: Optional[List[str]] = None) -> Material:
    """
    Center the material in the unit cell.

    Args:
        material (Material): The material object to center.
        axes (List[str]): The axes to center the material along.
    Returns:
        Material: The centered material object.
    """
    new_material = material.clone()
    new_material.to_crystal()
    if axes is None:
        axes = ["x", "y", "z"]
    min_x = get_atomic_coordinates_extremum(material, axis="x", extremum="min") if "x" in axes else 0
    max_x = get_atomic_coordinates_extremum(material, axis="x", extremum="max") if "x" in axes else 1
    min_y = get_atomic_coordinates_extremum(material, axis="y", extremum="min") if "y" in axes else 0
    max_y = get_atomic_coordinates_extremum(material, axis="y", extremum="max") if "y" in axes else 1
    if "z" in axes:
        material = translate_to_z_level(material, z_level="center")

    material = translate_by_vector(material, vector=[(1 - min_x - max_x) / 2, (1 - min_y - max_y) / 2, 0])
    return material


def wrap_to_unit_cell(material: Material) -> Material:
    """
    Wrap the material to the unit cell.

    Args:
        material (Material): The material to wrap.
    Returns:
        Material: The wrapped material.
    """
    atoms = to_ase(material)
    atoms.wrap()
    return Material(from_ase(atoms))


def filter_material_by_ids(material: Material, ids: List[int], invert: bool = False) -> Material:
    """
    Filter out only atoms corresponding to the ids.

    Args:
        material (Material): The material object to filter.
        ids (List[int]): The ids to filter by.
        invert (bool): Whether to invert the selection.

    Returns:
        Material: The filtered material object.
    """
    new_material = material.clone()
    new_basis = new_material.basis
    if invert is True:
        ids = list(set(new_basis.elements.ids) - set(ids))
    new_basis.filter_atoms_by_ids(ids)
    new_material.basis = new_basis
    return new_material


def filter_by_condition_on_coordinates(
    material: Material,
    condition: Callable[[List[float]], bool],
    use_cartesian_coordinates: bool = False,
    invert_selection: bool = False,
) -> Material:
    """
    Filter atoms based on a condition on their coordinates.

    Args:
        material (Material): The material object to filter.
        condition (Callable): The condition on coordinate function.
        use_cartesian_coordinates (bool): Whether to use cartesian coordinates.
        invert_selection (bool): Whether to invert the selection.

    Returns:
        Material: The filtered material object.
    """
    new_material = material.clone()
    ids = get_atom_indices_with_condition_on_coordinates(
        material,
        condition,
        use_cartesian_coordinates=use_cartesian_coordinates,
    )

    new_material = filter_material_by_ids(new_material, ids, invert=invert_selection)
    return new_material


def filter_by_layers(
    material: Material,
    center_coordinate: List[float] = [0, 0, 0],
    central_atom_id: Optional[int] = None,
    layer_thickness: float = 1.0,
    invert_selection: bool = False,
) -> Material:
    """
    Filter out atoms within a specified layer thickness of a central atom along c-vector direction.

    Args:
        material (Material): The material object to filter.
        center_coordinate (List[float]): Index of the central atom.
        central_atom_id (int): Index of the central atom.
        layer_thickness (float): Thickness of the layer in angstroms.
        invert_selection (bool): Whether to invert the selection.

    Returns:
        Material: The filtered material object.
    """
    if central_atom_id is not None:
        center_coordinate = material.basis.coordinates.get_element_value_by_index(central_atom_id)
    vectors = material.lattice.vectors
    direction_vector = vectors[2]

    def condition(coordinate):
        return is_coordinate_within_layer(coordinate, center_coordinate, direction_vector, layer_thickness)

    return filter_by_condition_on_coordinates(material, condition, invert_selection=invert_selection)


def filter_by_sphere(
    material: Material,
    center_coordinate: List[float] = [0, 0, 0],
    central_atom_id: Optional[int] = None,
    radius: float = 1,
    invert: bool = False,
) -> Material:
    """
    Filter out atoms within a specified radius of a central atom considering periodic boundary conditions.

    Args:
        material (Material): The material object to filter.
        central_atom_id (int): Index of the central atom.
        radius (float): Radius of the sphere in angstroms.
        invert (bool): Whether to invert the selection.

    Returns:
        Material: The filtered material object.
    """
    ids = get_atom_indices_within_radius_pbc(
        material=material,
        atom_index=central_atom_id,
        coordinate=center_coordinate,
        radius=radius,
    )
    return filter_material_by_ids(material, ids, invert=invert)


def filter_by_circle_projection(
    material: Material,
    x: float = 0.5,
    y: float = 0.5,
    r: float = 0.25,
    use_cartesian_coordinates: bool = False,
    invert_selection: bool = False,
) -> Material:
    """
    Get material with atoms that are within or outside an XY circle projection.

    Args:
        material (Material): The material object to filter.
        x (float): The x-coordinate of the circle center.
        y (float): The y-coordinate of the circle center.
        r (float): The radius of the circle.
        use_cartesian_coordinates (bool): Whether to use cartesian coordinates
        invert_selection (bool): Whether to invert the selection.

    Returns:
        Material: The filtered material object.
    """

    def condition(coordinate):
        return is_coordinate_in_cylinder(coordinate, [x, y, 0], r, min_z=0, max_z=1)

    return filter_by_condition_on_coordinates(
        material, condition, use_cartesian_coordinates=use_cartesian_coordinates, invert_selection=invert_selection
    )


def filter_by_cylinder(
    material: Material,
    center_position: List[float] = [0.5, 0.5],
    min_z: float = 0,
    max_z: float = 1,
    radius: float = 0.25,
    use_cartesian_coordinates: bool = False,
    invert_selection: bool = False,
) -> Material:
    """
    Get material with atoms that are within or outside a cylinder.

    Args:
        material (Material): The material object to filter.
        center_position (List[float]): The coordinates of the center position.
        radius (float): The radius of the cylinder.
        min_z (float): Lower limit of z-coordinate.
        max_z (float): Upper limit of z-coordinate.
        use_cartesian_coordinates (bool): Whether to use cartesian coordinates
        invert_selection (bool): Whether to invert the selection.

    Returns:
        Material: The filtered material object.
    """

    def condition(coordinate):
        return is_coordinate_in_cylinder(coordinate, center_position, radius, min_z, max_z)

    return filter_by_condition_on_coordinates(
        material, condition, use_cartesian_coordinates=use_cartesian_coordinates, invert_selection=invert_selection
    )


def filter_by_rectangle_projection(
    material: Material,
    min_coordinate: List[float] = [0, 0],
    max_coordinate: List[float] = [1, 1],
    use_cartesian_coordinates: bool = False,
    invert_selection: bool = False,
) -> Material:
    """
    Get material with atoms that are within or outside an XY rectangle projection.

    Args:
        material (Material): The material object to filter.
        min_coordinate (List[float]): The minimum coordinate of the rectangle.
        max_coordinate (List[float]): The maximum coordinate of the rectangle.
        use_cartesian_coordinates (bool): Whether to use cartesian coordinates
        invert_selection (bool): Whether to invert the selection.

    Returns:
        Material: The filtered material object.
    """
    min_coordinate = min_coordinate[:2] + [0]
    max_coordinate = max_coordinate[:2] + [1]

    def condition(coordinate):
        return is_coordinate_in_box(coordinate, min_coordinate, max_coordinate)

    return filter_by_condition_on_coordinates(
        material, condition, use_cartesian_coordinates=use_cartesian_coordinates, invert_selection=invert_selection
    )


def filter_by_box(
    material: Material,
    min_coordinate: List[float] = [0.0, 0.0, 0.0],
    max_coordinate: List[float] = [1.0, 1.0, 1.0],
    use_cartesian_coordinates: bool = False,
    invert_selection: bool = False,
) -> Material:
    """
    Get material with atoms that are within or outside an XYZ box.
    """

    def condition(coordinate):
        return is_coordinate_in_box(coordinate, min_coordinate, max_coordinate)

    return filter_by_condition_on_coordinates(
        material, condition, use_cartesian_coordinates=use_cartesian_coordinates, invert_selection=invert_selection
    )


def filter_by_triangle_projection(
    material: Material,
    coordinate_1: List[float] = [0, 0],
    coordinate_2: List[float] = [0, 1],
    coordinate_3: List[float] = [1, 0],
    min_z: float = 0,
    max_z: float = 1,
    use_cartesian_coordinates: bool = False,
    invert_selection: bool = False,
) -> Material:
    """
    Get material with atoms that are within or outside a prism formed by triangle projection.

    Args:
        material (Material): The material object to filter.
        coordinate_1 (List[float]): The coordinate of the first vertex.
        coordinate_2 (List[float]): The coordinate of the second vertex.
        coordinate_3 (List[float]): The coordinate of the third vertex.
        min_z (float): Lower limit of z-coordinate.
        max_z (float): Upper limit of z-coordinate.
        use_cartesian_coordinates (bool): Whether to use cartesian coordinates
        invert_selection (bool): Whether to invert the selection.

    Returns:
        Material: The filtered material object.
    """

    def condition(coordinate):
        return is_coordinate_in_triangular_prism(coordinate, coordinate_1, coordinate_2, coordinate_3, min_z, max_z)

    return filter_by_condition_on_coordinates(
        material, condition, use_cartesian_coordinates=use_cartesian_coordinates, invert_selection=invert_selection
    )


def add_vacuum(material: Material, vacuum: float = 5.0, on_top=True, to_bottom=False) -> Material:
    """
    Add vacuum to the material along the c-axis.
    On top, on bottom, or both.

    Args:
        material (Material): The material object to add vacuum to.
        vacuum (float): The thickness of the vacuum to add in angstroms.
        on_top (bool): Whether to add vacuum on top.
        to_bottom (bool): Whether to add vacuum on bottom.

    Returns:
        Material: The material object with vacuum added.
    """
    new_material_atoms = to_ase(material)
    vacuum_amount = vacuum * 2 if on_top and to_bottom else vacuum
    ase_add_vacuum(new_material_atoms, vacuum_amount)
    new_material = Material(from_ase(new_material_atoms))
    if to_bottom and not on_top:
        new_material = translate_to_z_level(new_material, z_level="top")
    elif on_top and to_bottom:
        new_material = translate_to_z_level(new_material, z_level="center")
    return new_material


def remove_vacuum(material: Material, from_top=True, from_bottom=True, fixed_padding=1.0) -> Material:
    """
    Remove vacuum from the material along the c-axis.
    From top, from bottom, or from both.

    Args:
        material (Material): The material object to set the vacuum thickness.
        from_top (bool): Whether to remove vacuum from the top.
        from_bottom (bool): Whether to remove vacuum from the bottom.
        fixed_padding (float): The fixed padding of vacuum to add to avoid collisions in pbc (in angstroms).

    Returns:
        Material: The material object with the vacuum thickness set.
    """
    translated_material = translate_to_z_level(material, z_level="bottom")
    new_basis = translated_material.basis
    new_basis.to_cartesian()
    new_lattice = translated_material.lattice
    new_lattice.c = get_atomic_coordinates_extremum(translated_material, use_cartesian_coordinates=True) + fixed_padding
    new_basis.cell.vector3 = new_lattice.vectors[2]
    new_basis.to_crystal()
    new_material = material.clone()

    new_material.basis = new_basis
    new_material.lattice = new_lattice

    if from_top and not from_bottom:
        new_material = translate_to_z_level(new_material, z_level="top")
    if from_bottom and not from_top:
        new_material = translate_to_z_level(new_material, z_level="bottom")
    return new_material


def rotate_material(material: Material, axis: List[int], angle: float) -> Material:
    """
    Rotate the material around a given axis by a specified angle.

    Args:
        material (Material): The material to rotate.
        axis (List[int]): The axis to rotate around, expressed as [x, y, z].
        angle (float): The angle of rotation in degrees.
    Returns:
        Atoms: The rotated material.
    """
    if material.basis.is_in_cartesian_units:
        crystal_basis = material.basis.copy()
        crystal_basis.to_crystal()
        material.basis = crystal_basis
    atoms = to_ase(material)
    atoms.rotate(v=axis, a=angle, center="COM")
    atoms.wrap()

    return Material(from_ase(atoms))


def get_passivant_coordinate_crystal(
    material: Material, atom_coordinate: np.ndarray, bond_direction: np.ndarray, bond_length: float = 2.0
):
    bond_normal_cartesian = np.array(material.basis.cell.convert_point_to_cartesian(bond_direction.tolist()))
    bond_normal_cartesian /= np.linalg.norm(bond_normal_cartesian)
    bond_vector_cartesian = bond_normal_cartesian * bond_length
    bond_vector_crystal = np.array(material.basis.cell.convert_point_to_crystal(bond_vector_cartesian.tolist()))
    passivant_coordinate_crystal = atom_coordinate + bond_vector_crystal
    return passivant_coordinate_crystal


def passivate_material(
    slab: Material, passivant: str, bond_length: float = 2.0, coordination_threshold: Optional[int] = None
):
    nudge_value = 0.01
    supercell_scaling_factor = [3, 3, 3]
    min_coordinate = [1 / 3, 1 / 3, 1 / 3]
    max_coordinate = [2 / 3, 2 / 3, 2 / 3]
    adjusted_min_coordinate = (np.array(min_coordinate) - nudge_value).tolist()
    adjusted_max_coordinate = (np.array(max_coordinate) + nudge_value).tolist()
    slab = translate_to_z_level(slab, "center")
    new_basis = slab.basis.copy()
    slab_supercell = create_supercell(slab, scaling_factor=supercell_scaling_factor)

    central_cell_ids = [
        id
        for id, coordinate in zip(slab_supercell.basis.coordinates.ids, slab_supercell.basis.coordinates.values)
        if is_coordinate_in_box(coordinate, adjusted_min_coordinate, adjusted_max_coordinate)
    ]
    undercoordinated_atom_indices, atom_neighbors_info = get_undercoordinated_atoms(
        slab_supercell, central_cell_ids, coordination_threshold
    )
    for index in undercoordinated_atom_indices:
        atom_coordinate_crystal = np.array(slab_supercell.basis.coordinates.values[index])
        neighbors_average_coordinate_crystal = atom_neighbors_info[index][2]
        bond_normal_crystal = np.array(slab_supercell.basis.coordinates.values[index]) - np.array(
            neighbors_average_coordinate_crystal
        )
        passivant_coordinate_crystal = get_passivant_coordinate_crystal(
            slab_supercell, atom_coordinate_crystal, bond_normal_crystal, bond_length
        ).tolist()
        passivant_coordinate_crystal_original_cell = transform_coordinate_to_supercell(
            passivant_coordinate_crystal,
            scaling_factor=supercell_scaling_factor,
            translation_vector=min_coordinate,
            reverse=True,
        )
        new_basis.add_atom(passivant, passivant_coordinate_crystal_original_cell)

    slab.basis = new_basis
    return slab


# TODO: Get this from peroidic table
BOND_LENGTHS_MAP = {
    ("C", "H"): 1.09,
    ("Ni", "H"): 1.09,
    ("Si", "H"): 1.48,
}


class SURFACE_TYPES(str, Enum):
    TOP = "top"
    BOTTOM = "bottom"
    BOTH = "both"


def add_passivant_atoms_to_material(
    material: Material,
    axis: Literal["x", "y", "z"],
    surface: Optional[SURFACE_TYPES] = SURFACE_TYPES.BOTH,
    tolerance: float = 0.01,
    passivant: Optional[str] = "H",
    default_bond_length: float = 1.0,
) -> Material:
    """
    Add passivant atoms to the specified surface or edge of the material.

    Args:
        material (Material): The material object to add passivant atoms to.
        axis (AXIS_TYPES): The axis along which the surface or edge lies ("x", "y", or "z").
        surface (SURFACE_TYPES): The surface to add passivant atoms to ("top", "bottom", or "both"). For edges, assumes both.
        tolerance (float): The tolerance for selecting surface atoms.
        passivant (str): The chemical symbol of the passivating atom (e.g., 'H').
        default_bond_length (float): The default bond length to use if the pair is not found in BOND_LENGTHS_MAP.

    Returns:
        Material: The material object with passivation atoms added.
    """
    axis_idx = {"x": 0, "y": 1, "z": 2}[axis]
    basis = material.basis.copy()
    if axis == "z":
        if surface in [SURFACE_TYPES.TOP, SURFACE_TYPES.BOTH]:
            z_extremum = get_atomic_coordinates_extremum(material, "max", axis)
            add_passivant_atoms_to_basis(
                material, basis, axis_idx, z_extremum, tolerance, passivant, default_bond_length, positive=True
            )
        if surface in [SURFACE_TYPES.BOTTOM, SURFACE_TYPES.BOTH]:
            z_extremum = get_atomic_coordinates_extremum(material, "min", axis)
            add_passivant_atoms_to_basis(
                material, basis, axis_idx, z_extremum, tolerance, passivant, default_bond_length, positive=False
            )
    else:
        x_or_y_extremum_max = get_atomic_coordinates_extremum(material, "max", axis)
        add_passivant_atoms_to_basis(
            material, basis, axis_idx, x_or_y_extremum_max, tolerance, passivant, default_bond_length, positive=True
        )

        x_or_y_extremum_min = get_atomic_coordinates_extremum(material, "min", axis)
        add_passivant_atoms_to_basis(
            material, basis, axis_idx, x_or_y_extremum_min, tolerance, passivant, default_bond_length, positive=False
        )
    material.basis = basis
    return material


def add_passivant_atoms_to_basis(
    material, basis, axis_idx, extremum, tolerance, passivant, default_bond_length, positive=True
):
    """
    Helper function to add passivant atoms to the specified extremum along a given axis.
    """
    min_box = [0, 0, 0]
    max_box = [1, 1, 1]
    min_box[axis_idx] = extremum - tolerance
    max_box[axis_idx] = extremum + tolerance

    surface_atoms = filter_by_box(material, min_box, max_box)
    surface_indices = surface_atoms.basis.coordinates.ids

    for idx in surface_indices:
        atom_coordinate = material.basis.coordinates.values[idx]
        element = material.basis.elements.values[idx]
        bond_length = BOND_LENGTHS_MAP.get((element, passivant), default_bond_length)

        bond_length_vector = [0, 0, 0]
        bond_length_vector[axis_idx] = bond_length if positive else -bond_length
        bond_length_crystal = material.basis.cell.convert_point_to_crystal(bond_length_vector)

        passivant_coordinate = atom_coordinate + bond_length_crystal
        basis.add_atom(passivant, passivant_coordinate)


def passivate_surface(
    material: Material,
    passivant: str = "H",
    default_bond_length: float = 1.0,
    surface: SURFACE_TYPES = SURFACE_TYPES.BOTH,
) -> Material:
    """
    Passivates the top and/or bottom surfaces of a material by adding atoms along the Z-axis,
    with bond lengths determined by the element and passivant.

    Args:
        material (Material): The material to passivate.
        passivant (str): The chemical symbol of the passivating atom (e.g., 'H').
        default_bond_length (float): The default bond length to use if the pair is not found in BOND_LENGTHS_MAP.
        surface (SURFACE_TYPES): The surface to passivate ("top", "bottom", or "both").

    Returns:
        Material: The passivated material.
    """
    material = translate_to_z_level(material, "center")

    return add_passivant_atoms_to_material(
        material=material,
        axis="z",
        surface=surface,
        passivant=passivant,
        default_bond_length=default_bond_length,
    )


def passivate_edges(
    material: Material,
    passivant: str = "H",
    default_bond_length: float = 1.0,
    axis: Literal["x", "y"] = "x",
) -> Material:
    """
    Passivates the edges of a 2D material by adding atoms along the X or Y axis,
    with bond lengths determined by the element and passivant.

    Args:
        material (Material): The material to passivate.
        passivant (str): The chemical symbol of the passivating atom (e.g., 'H').
        default_bond_length (float): The default bond length to use if the pair is not found in BOND_LENGTHS_MAP.
        axis (AXIS_TYPES): The axis along which the edges lie ("x" or "y").

    Returns:
        Material: The passivated material.
    """
    material = translate_to_z_level(material, "center")

    return add_passivant_atoms_to_material(
        material=material, axis=axis, passivant=passivant, default_bond_length=default_bond_length
    )
