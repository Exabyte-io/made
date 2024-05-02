from ase import Atoms
from ase.calculators.calculator import Calculator

from .analyze import get_surface_area
from .convert import decorator_convert_material_args_kwargs_to_atoms


@decorator_convert_material_args_kwargs_to_atoms
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


@decorator_convert_material_args_kwargs_to_atoms
def calculate_total_energy_per_atom(atoms: Atoms, calculator: Calculator):
    """
    Set calculator for ASE Atoms and calculate the total energy per atom.

    Args:
            atoms (ase.Atoms): The Atoms object to calculate the energy of.
            calculator (ase.calculators.calculator.Calculator): The calculator to use for the energy calculation.

    Returns:
            float: The energy per atom of the atoms.
    """
    return calculate_total_energy(atoms, calculator) / atoms.get_global_number_of_atoms()


@decorator_convert_material_args_kwargs_to_atoms
def calculate_surface_energy(slab: Atoms, bulk: Atoms, calculator: Calculator):
    """
    Calculate the surface energy by subtracting the weighted bulk energy from the slab energy.

    Args:
        slab (ase.Atoms): The slab Atoms object to calculate the surface energy of.
        bulk (ase.Atoms): The bulk Atoms object to calculate the surface energy of.
        calculator (ase.calculators.calculator.Calculator): The calculator to use for the energy calculation.

    Returns:
        float: The surface energy of the slab.
    """
    number_of_atoms = slab.get_global_number_of_atoms()
    area = get_surface_area(slab)
    return (
        calculate_total_energy(slab, calculator) - calculate_total_energy_per_atom(bulk, calculator) * number_of_atoms
    ) / (2 * area)


@decorator_convert_material_args_kwargs_to_atoms
def calculate_adhesion_energy(interface: Atoms, substrate_slab: Atoms, layer_slab: Atoms, calculator: Calculator):
    """
    Calculate the adhesion energy.
    The adhesion energy is the difference between the energy of the interface and
    the sum of the energies of the substrate and layer.
    According to: 10.1088/0953-8984/27/30/305004

    Args:
        interface (ase.Atoms): The interface Atoms object to calculate the adhesion energy of.
        substrate_slab (ase.Atoms): The substrate slab Atoms object to calculate the adhesion energy of.
        layer_slab (ase.Atoms): The layer slab Atoms object to calculate the adhesion energy of.
        calculator (ase.calculators.calculator.Calculator): The calculator to use for the energy calculation.

    Returns:
        float: The adhesion energy of the interface.
    """
    energy_substrate_slab = calculate_total_energy(substrate_slab, calculator)
    energy_layer_slab = calculate_total_energy(layer_slab, calculator)
    energy_interface = calculate_total_energy(interface, calculator)
    area = get_surface_area(interface)
    return (energy_substrate_slab + energy_layer_slab - energy_interface) / area


@decorator_convert_material_args_kwargs_to_atoms
def calculate_interfacial_energy(
    interface: Atoms,
    substrate_slab: Atoms,
    substrate_bulk: Atoms,
    layer_slab: Atoms,
    layer_bulk: Atoms,
    calculator: Calculator,
):
    """
    Calculate the interfacial energy.
    The interfacial energy is the sum of the surface energies of the substrate and layer minus the adhesion energy.
    According to Dupré's formula

    Args:
        interface (ase.Atoms): The interface Atoms object to calculate the interfacial energy of.
        substrate_slab (ase.Atoms): The substrate slab Atoms object to calculate the interfacial energy of.
        substrate_bulk (ase.Atoms): The substrate bulk Atoms object to calculate the interfacial energy of.
        layer_slab (ase.Atoms): The layer slab Atoms object to calculate the interfacial energy of.
        layer_bulk (ase.Atoms): The layer bulk Atoms object to calculate the interfacial energy of.
        calculator (ase.calculators.calculator.Calculator): The calculator to use for the energy calculation.

    Returns:
        float: The interfacial energy of the interface.
    """

    surface_energy_substrate = calculate_surface_energy(substrate_slab, substrate_bulk, calculator)
    surface_energy_layer = calculate_surface_energy(layer_slab, layer_bulk, calculator)
    adhesion_energy = calculate_adhesion_energy(interface, substrate_slab, layer_slab, calculator)
    return surface_energy_layer + surface_energy_substrate - adhesion_energy
