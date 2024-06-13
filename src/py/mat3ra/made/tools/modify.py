from typing import Union, List

import numpy as np
from mat3ra.made.material import Material
from mat3ra.made.utils import filter_array_with_id_value_by_ids, filter_array_with_id_value_by_values
from pymatgen.analysis.structure_analyzer import SpacegroupAnalyzer
from pymatgen.core.structure import Structure

from .analyze import select_atoms_within_layers, select_atoms_within_radius_pbc
from .convert import decorator_convert_material_args_kwargs_to_structure
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
    labels = material.basis["labels"]
    filtered_labels = filter_array_with_id_value_by_values(labels, label)
    filtered_label_ids = [item["id"] for item in filtered_labels]
    for key in ["coordinates", "elements", "labels"]:
        new_material.basis[key] = filter_array_with_id_value_by_ids(new_material.basis[key], filtered_label_ids)
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


def filter_by_layers(material, central_atom_id, layer_thickness, invert=False):
    new_material = material.clone()
    ids = select_atoms_within_layers(
        material,
        central_atom_id,
        layer_thickness,
    )
    if invert:
        ids = [i for i in range(len(material.basis["coordinates"])) if i not in ids]
    new_basis = filter_array_with_id_value_by_ids(material.basis, ids)
    new_material.basis = new_basis
    return new_material


def filter_by_sphere(material, central_atom_id, radius):
    new_material = material.clone()
    ids = select_atoms_within_radius_pbc(
        material,
        central_atom_id,
        radius,
    )
    return new_material
