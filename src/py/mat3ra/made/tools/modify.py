from enum import Enum
from typing import Union

from ase import Atoms
from ase.optimize import BFGS
from pymatgen.analysis.structure_analyzer import SpacegroupAnalyzer
from pymatgen.core.structure import Structure

from .calculate import CalculatorEnum, calculator_by_name_map
from .convert import (
    decorator_convert_material_args_kwargs_to_atoms,
    decorator_convert_material_args_kwargs_to_structure,
)
from .utils import translate_to_bottom_pymatgen_structure


# TODO: ASE related enums and maps should be placed close to each other
class OptimizerEnum(Enum):
    BFGS = "BFGS"


optimizer_by_name_map = {"BFGS": BFGS}


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


@decorator_convert_material_args_kwargs_to_atoms
def relax_atoms(atoms: Atoms, calculator: CalculatorEnum, optimizer: OptimizerEnum, **kwargs):
    """
    Relax the atoms using the calculator.

    Args:
        atoms (ase.Atoms): The Atoms object to relax.
        calculator (CalculatorEnum): The calculator name to use for the relaxation.
        optimizer (OptimizerEnum): The optimizer name to use for the relaxation.
        **kwargs: Additional keyword arguments to pass to the relaxation method.

    Returns:
        ase.Atoms: The relaxed Atoms object.
    """
    calculator_object = calculator_by_name_map[calculator.value]()
    atoms.set_calculator(calculator_object)
    atoms.get_potential_energy()

    optimizer_object = optimizer_by_name_map[optimizer.value]
    dyn = optimizer_object(atoms, **kwargs)
    dyn.run(fmax=0.01)
    return atoms
