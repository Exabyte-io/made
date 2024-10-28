from typing import Callable, List, Literal, Optional, Tuple, Union

import numpy as np
from mat3ra.made.material import Material

from .analyze import (
    get_atom_indices_with_condition_on_coordinates,
    get_atom_indices_within_radius_pbc,
    get_atomic_coordinates_extremum,
)
from .convert import from_ase, to_ase
from .convert.utils import InterfacePartsEnum
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


def filter_by_ids(material: Material, ids: List[int], invert: bool = False) -> Material:
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

    new_material = filter_by_ids(new_material, ids, invert=invert_selection)
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


def get_default_min_max(material: Material, use_cartesian_coordinates: bool) -> Tuple[List[float], List[float]]:
    """
    Get default min and max coordinates for the material based on the coordinate system (crystal or cartesian).

    Args:
        material (Material): The material object.
        use_cartesian_coordinates (bool): Whether to use Cartesian coordinates.

    Returns:
        Tuple[List[float], List[float]]: Default min and max coordinates.
    """
    axes: List[Literal["x", "y", "z"]] = ["x", "y", "z"]
    min_coords = []
    max_coords = []
    for axis in axes:
        min_val = get_atomic_coordinates_extremum(
            material, "min", axis, use_cartesian_coordinates=use_cartesian_coordinates
        )
        max_val = get_atomic_coordinates_extremum(
            material, "max", axis, use_cartesian_coordinates=use_cartesian_coordinates
        )
        min_coords.append(min_val)
        max_coords.append(max_val)
    return min_coords, max_coords


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
    return filter_by_ids(material, ids, invert=invert)


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
    center_position: Optional[List[float]] = None,
    min_z: Optional[float] = None,
    max_z: Optional[float] = None,
    radius: float = 0.25,
    use_cartesian_coordinates: bool = False,
    invert_selection: bool = False,
) -> Material:
    """
    Get material with atoms that are within or outside a cylinder.

    Args:
        material (Material): The material object to filter.
        center_position (List[float], optional): The center position of the cylinder. Defaults to material's center.
        min_z (float, optional): Lower limit of z-coordinate. Defaults to material's min z.
        max_z (float, optional): Upper limit of z-coordinate. Defaults to material's max z.
        radius (float): The radius of the cylinder.
        use_cartesian_coordinates (bool): Whether to use cartesian coordinates.
        invert_selection (bool): Whether to invert the selection.

    Returns:
        Material: The filtered material object.
    """
    if center_position is None:
        default_min, default_max = get_default_min_max(material, use_cartesian_coordinates)
        center_position = [(min_c + max_c) / 2 for min_c, max_c in zip(default_min, default_max)]

    if min_z is None or max_z is None:
        min_z = get_atomic_coordinates_extremum(
            material, "min", "z", use_cartesian_coordinates=use_cartesian_coordinates
        )
        max_z = get_atomic_coordinates_extremum(
            material, "max", "z", use_cartesian_coordinates=use_cartesian_coordinates
        )

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
    min_z = get_atomic_coordinates_extremum(material, "min", "z", use_cartesian_coordinates=use_cartesian_coordinates)
    max_z = get_atomic_coordinates_extremum(material, "max", "z", use_cartesian_coordinates=use_cartesian_coordinates)
    min_coordinate = min_coordinate[:2] + [min_z]
    max_coordinate = max_coordinate[:2] + [max_z]

    def condition(coordinate):
        return is_coordinate_in_box(coordinate, min_coordinate, max_coordinate)

    return filter_by_condition_on_coordinates(
        material, condition, use_cartesian_coordinates=use_cartesian_coordinates, invert_selection=invert_selection
    )


def filter_by_box(
    material: Material,
    min_coordinate: Optional[List[float]] = None,
    max_coordinate: Optional[List[float]] = None,
    use_cartesian_coordinates: bool = False,
    invert_selection: bool = False,
) -> Material:
    """
    Get material with atoms that are within or outside an XYZ box.

    Args:
        material (Material): The material to filter.
        min_coordinate (List[float], optional): The minimum coordinate of the box. Defaults to material's min.
        max_coordinate (List[float], optional): The maximum coordinate of the box. Defaults to material's max.
        use_cartesian_coordinates (bool): Whether to use cartesian coordinates.
        invert_selection (bool): Whether to invert the selection.

    Returns:
        Material: The filtered material object.
    """
    if min_coordinate is None or max_coordinate is None:
        default_min, default_max = get_default_min_max(material, use_cartesian_coordinates)
        min_coordinate = min_coordinate if min_coordinate is not None else default_min
        max_coordinate = max_coordinate if max_coordinate is not None else default_max

    def condition(coordinate):
        return is_coordinate_in_box(coordinate, min_coordinate, max_coordinate)

    return filter_by_condition_on_coordinates(
        material, condition, use_cartesian_coordinates=use_cartesian_coordinates, invert_selection=invert_selection
    )


def filter_by_triangle_projection(
    material: Material,
    coordinate_1: Optional[List[float]] = None,
    coordinate_2: Optional[List[float]] = None,
    coordinate_3: Optional[List[float]] = None,
    min_z: Optional[float] = None,
    max_z: Optional[float] = None,
    use_cartesian_coordinates: bool = False,
    invert_selection: bool = False,
) -> Material:
    """
    Get material with atoms that are within or outside a triangular prism.

    Args:
        material (Material): The material object to filter.
        coordinate_1 (List[float], optional): First vertex of the triangle. Defaults to material's corner.
        coordinate_2 (List[float], optional): Second vertex of the triangle. Defaults to material's corner.
        coordinate_3 (List[float], optional): Third vertex of the triangle. Defaults to material's corner.
        min_z (float, optional): Lower z-limit. Defaults to material's min z.
        max_z (float, optional): Upper z-limit. Defaults to material's max z.
        use_cartesian_coordinates (bool): Whether to use cartesian coordinates.
        invert_selection (bool): Whether to invert the selection.

    Returns:
        Material: The filtered material object.
    """
    if coordinate_1 is None or coordinate_2 is None or coordinate_3 is None:
        default_min, default_max = get_default_min_max(material, use_cartesian_coordinates)
        coordinate_1 = coordinate_1 if coordinate_1 is not None else [default_min[0], default_min[1]]
        coordinate_2 = coordinate_2 if coordinate_2 is not None else [default_min[0], default_max[1]]
        coordinate_3 = coordinate_3 if coordinate_3 is not None else [default_max[0], default_min[1]]

    if min_z is None:
        min_z = get_atomic_coordinates_extremum(
            material, "min", "z", use_cartesian_coordinates=use_cartesian_coordinates
        )
    if max_z is None:
        max_z = get_atomic_coordinates_extremum(
            material, "max", "z", use_cartesian_coordinates=use_cartesian_coordinates
        )

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


def add_vacuum_sides(material: Material, vacuum: float = 5.0, on_x=False, on_y=False) -> Material:
    """
    Add vacuum to the material along the x-axis and/or y-axis.
    On x, on y, or both.

    Args:
        material (Material): The material object to add vacuum to.
        vacuum (float): The thickness of the vacuum to add in angstroms.
        on_x (bool): Whether to add vacuum on the x-axis.
        on_y (bool): Whether to add vacuum on the y-axis.

    Returns:
        Material: The material object with vacuum added.
    """
    new_material = material.clone()
    min_x = get_atomic_coordinates_extremum(new_material, "min", "x", use_cartesian_coordinates=True)
    max_x = get_atomic_coordinates_extremum(new_material, "max", "x", use_cartesian_coordinates=True)
    min_y = get_atomic_coordinates_extremum(new_material, "min", "y", use_cartesian_coordinates=True)
    max_y = get_atomic_coordinates_extremum(new_material, "max", "y", use_cartesian_coordinates=True)

    new_lattice_a_vector, new_lattice_b_vector, new_lattice_c_vector = new_material.lattice.vector_arrays
    if on_x:
        new_x_length = max_x - min_x + 2 * vacuum
        new_lattice_a_vector = [new_x_length, 0, 0]
    if on_y:
        new_y_length = max_y - min_y + 2 * vacuum
        new_lattice_b_vector = [0, new_y_length, 0]

    new_material.set_new_lattice_vectors(new_lattice_a_vector, new_lattice_b_vector, new_lattice_c_vector)
    new_material = translate_by_vector(
        new_material,
        [-min_x + vacuum if on_x else 0, -min_y + vacuum if on_y else 0, 0],
        use_cartesian_coordinates=True,
    )

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


def rotate(material: Material, axis: List[int], angle: float, wrap: bool = True, rotate_cell=False) -> Material:
    """
    Rotate the material around a given axis by a specified angle.

    Args:
        material (Material): The material to rotate.
        axis (List[int]): The axis to rotate around, expressed as [x, y, z].
        angle (float): The angle of rotation in degrees.
        wrap (bool): Whether to wrap the material to the unit cell.
        rotate_cell (bool): Whether to rotate the cell.
    Returns:
        Atoms: The rotated material.
    """
    if material.basis.is_in_cartesian_units:
        crystal_basis = material.basis.copy()
        crystal_basis.to_crystal()
        material.basis = crystal_basis
    atoms = to_ase(material)
    atoms.rotate(v=axis, a=angle, center="COU", rotate_cell=rotate_cell)
    if wrap:
        atoms.wrap()

    return Material(from_ase(atoms))


def interface_displace_part(
    interface: Material,
    displacement: List[float],
    label: InterfacePartsEnum = InterfacePartsEnum.FILM,
    use_cartesian_coordinates=True,
) -> Material:
    """
    Displace atoms in an interface along a certain direction.

    Args:
        interface (Material): The interface Material object.
        displacement (List[float]): The displacement vector in angstroms or crystal coordinates.
        label (InterfacePartsEnum): The label of the atoms to displace ("substrate" or "film").
        use_cartesian_coordinates (bool): Whether to use cartesian coordinates.

    Returns:
        Material: The displaced material object.
    """
    new_material = interface.clone()
    if use_cartesian_coordinates:
        new_material.to_cartesian()
    labels_array = new_material.basis.labels.to_array_of_values_with_ids()
    displaced_label_ids = [_label.id for _label in labels_array if _label.value == int(label)]

    new_coordinates_values = new_material.basis.coordinates.values
    for atom_id in displaced_label_ids:
        current_coordinate = new_material.basis.coordinates.get_element_value_by_index(atom_id)
        new_atom_coordinate = np.array(current_coordinate) + np.array(displacement)
        new_coordinates_values[atom_id] = new_atom_coordinate

    new_material.set_coordinates(new_coordinates_values)
    new_material.to_crystal()
    new_material = wrap_to_unit_cell(new_material)
    return new_material


def interface_get_part(
    interface: Material,
    part: InterfacePartsEnum = InterfacePartsEnum.FILM,
) -> Material:
    if interface.metadata["build"]["configuration"]["type"] != "InterfaceConfiguration":
        raise ValueError("The material is not an interface.")
    interface_part_material = interface.clone()
    film_atoms_basis = interface_part_material.basis.filter_atoms_by_labels([int(part)])
    interface_part_material.basis = film_atoms_basis
    return interface_part_material
