from typing import List, Union

from mat3ra.made.material import Material

from .analyze import get_atom_indices_within_layer_by_atom_index, get_atom_indices_within_radius_pbc
from .convert import decorator_convert_material_args_kwargs_to_structure
from .third_party import PymatgenSpacegroupAnalyzer, PymatgenStructure
from .utils import translate_to_bottom_pymatgen_structure


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
def translate_to_bottom(structure: PymatgenStructure, use_conventional_cell: bool = True):
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
        structure = PymatgenSpacegroupAnalyzer(structure).get_conventional_standard_structure()
    structure = translate_to_bottom_pymatgen_structure(structure)
    return structure


@decorator_convert_material_args_kwargs_to_structure
def wrap_to_unit_cell(structure: PymatgenStructure):
    """
    Wrap atoms to the cell

    Args:
        structure (PymatgenStructure): The pymatgen PymatgenStructure object to normalize.
    Returns:
        PymatgenStructure: The wrapped pymatgen PymatgenStructure object.
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
