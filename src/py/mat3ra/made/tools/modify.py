from typing import Union

from ase import Atoms

from .convert import decorator_convert_material_args_kwargs_to_atoms


@decorator_convert_material_args_kwargs_to_atoms
def filter_by_label(atoms: Atoms, label: Union[int, str]):
    """
    Filter out only atoms corresponding to the label/tag.

    Args:
        atoms (ase.Atoms): The Atoms object to filter.
        label (int|str): The tag/label to filter by.

    Returns:
        ase.Atoms: The filtered Atoms object.
    """
    return atoms[atoms.get_tags() == label]
