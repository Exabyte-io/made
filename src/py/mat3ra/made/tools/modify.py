from typing import Callable, List, Union

from mat3ra.made.material import Material
from pymatgen.analysis.structure_analyzer import SpacegroupAnalyzer
from pymatgen.core.structure import Structure

from .analyze import (
    get_atom_indices_with_coordinates_condition,
    get_atom_indices_within_layer_by_atom_index,
    get_atom_indices_within_radius_pbc,
)
from .convert import decorator_convert_material_args_kwargs_to_structure
from .utils import is_point_in_circle, is_point_in_rectangle, translate_to_bottom_pymatgen_structure


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


@decorator_convert_material_args_kwargs_to_structure
def translate_to_bottom(structure: Structure, use_conventional_cell: bool = True):
    """
    Translate atoms to the bottom of the cell (vacuum on top) to allow for the correct consecutive interface generation.
    If use_conventional_cell is passed, conventional cell is used.

    Args:
        structure (Structure): The pymatgen Structure object to normalize.
        use_conventional_cell: Whether to convert to the conventional cell.
    Returns:
        Structure: The normalized pymatgen Structure object.
    """
    if use_conventional_cell:
        structure = SpacegroupAnalyzer(structure).get_conventional_standard_structure()
    structure = translate_to_bottom_pymatgen_structure(structure)
    return structure


@decorator_convert_material_args_kwargs_to_structure
def wrap_to_unit_cell(structure: Structure):
    """
    Wrap atoms to the cell

    Args:
        structure (Structure): The pymatgen Structure object to normalize.
    Returns:
        Structure: The wrapped pymatgen Structure object.
    """
    structure.make_supercell((1, 1, 1), to_unit_cell=True)
    return structure


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


def filter_by_coordinates_condition(
    material: Material,
    condition: Callable[[List[float]], bool],
    use_cartesian_coordinates: bool = False,
    invert_selection: bool = False,
) -> Material:
    """
    Filter atoms based on a condition on their coordinates.

    Args:
        material (Material): The material object to filter.
        condition (Callable): The coordinate condition function.
        use_cartesian_coordinates (bool): Whether to use cartesian coordinates.
        invert_selection (bool): Whether to invert the selection.

    Returns:
        Material: The filtered material object.
    """
    new_material = material.clone()
    ids = get_atom_indices_with_coordinates_condition(
        material,
        condition,
        use_cartesian_coordinates=use_cartesian_coordinates,
    )

    new_material = filter_material_by_ids(new_material, ids, invert=invert_selection)
    return new_material


def filter_by_layers(
    material: Material, central_atom_id: int, layer_thickness: float, invert: bool = False
) -> Material:
    """
    Filter out atoms within a specified layer thickness of a central atom along c-vector direction.

    Args:
        material (Material): The material object to filter.
        central_atom_id (int): Index of the central atom.
        layer_thickness (float): Thickness of the layer in angstroms.
        invert (bool): Whether to invert the selection.

    Returns:
        Material: The filtered material object.
    """
    ids = get_atom_indices_within_layer_by_atom_index(
        material,
        central_atom_id,
        layer_thickness,
    )
    return filter_material_by_ids(material, ids, invert=invert)


def filter_by_sphere(material: Material, central_atom_id: int, radius: float, invert: bool = False) -> Material:
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
    condition = is_point_in_circle(x, y, r)
    return filter_by_coordinates_condition(
        material, condition, use_cartesian_coordinates=use_cartesian_coordinates, invert_selection=invert_selection
    )


def filter_by_rectangle_projection(
    material: Material,
    x_min: float = 0.25,
    x_max: float = 0.75,
    y_min: float = 0.25,
    y_max: float = 0.75,
    use_cartesian_coordinates: bool = False,
    invert_selection: bool = False,
) -> Material:
    """
    Get material with atoms that are within or outside an XY rectangle projection.

    Args:

        material (Material): The material object to filter.
        x_min (float): The minimum x-coordinate of the rectangle.
        x_max (float): The maximum x-coordinate of the rectangle.
        y_min (float): The minimum y-coordinate of the rectangle.
        y_max (float): The maximum y-coordinate of the rectangle.
        use_cartesian_coordinates (bool): Whether to use cartesian coordinates
        invert_selection (bool): Whether to invert the selection.

    Returns:
        Material: The filtered material object.
    """
    condition = is_point_in_rectangle(x_min, x_max, y_min, y_max)

    return filter_by_coordinates_condition(
        material, condition, use_cartesian_coordinates=use_cartesian_coordinates, invert_selection=invert_selection
    )
