from ase import Atoms
from ase.calculators.calculator import Calculator

from .convert import convert_material_args_kwargs_to_atoms


@convert_material_args_kwargs_to_atoms
def calculate_total_energy(atoms: Atoms, calculator: Calculator):
    """
    Set calculator for ASE Atoms and calculate the total energy.

    Args:
        atoms (ase.Atoms): The Atoms object to calculate the energy of.
        calculator (ase.calculators.calculator.Calculator): The calculator to use for the energy calculation.

    Returns:
        float: The total energy of the atoms.
    """
    atoms.set_calculator(calculator)
    return atoms.get_total_energy()
