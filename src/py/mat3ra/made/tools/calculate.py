from typing import Optional

from ase import Atoms
from ase.calculators.calculator import Calculator
from ase.calculators.emt import EMT

from ..material import Material
from .analyze import get_surface_area
from .build.interface.utils import get_slab
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
def calculate_adhesion_energy(interface: Atoms, substrate_slab: Atoms, film_slab: Atoms, calculator: Calculator):
    """
    Calculate the adhesion energy.
    The adhesion energy is the difference between the energy of the interface and
    the sum of the energies of the substrate and film.
    According to: 10.1088/0953-8984/27/30/305004

    Args:
        interface (ase.Atoms): The interface Atoms object to calculate the adhesion energy of.
        substrate_slab (ase.Atoms): The substrate slab Atoms object to calculate the adhesion energy of.
        film_slab (ase.Atoms): The film slab Atoms object to calculate the adhesion energy of.
        calculator (ase.calculators.calculator.Calculator): The calculator to use for the energy calculation.

    Returns:
        float: The adhesion energy of the interface.
    """
    energy_substrate_slab = calculate_total_energy(substrate_slab, calculator)
    energy_film_slab = calculate_total_energy(film_slab, calculator)
    energy_interface = calculate_total_energy(interface, calculator)
    area = get_surface_area(interface)
    return (energy_substrate_slab + energy_film_slab - energy_interface) / area


def calculate_interfacial_energy(
    interface: Material,
    substrate_slab: Optional[Material] = None,
    substrate_bulk: Optional[Material] = None,
    film_slab: Optional[Material] = None,
    film_bulk: Optional[Material] = None,
    calculator: Calculator = EMT(),
):
    """
    Calculate the interfacial energy.
    The interfacial energy is the sum of the surface energies of the substrate and film minus the adhesion energy.
    According to Dupré's formula

    Args:
        interface (Material): The interface Material object to calculate the interfacial energy of.
        substrate_slab (Material): The substrate slab Material object to calculate the interfacial energy of.
        substrate_bulk (Material): The substrate bulk Material object to calculate the interfacial energy of.
        film_slab (Material): The film slab Material object to calculate the interfacial energy of.
        film_bulk (Material): The film bulk Material object to calculate the interfacial energy of.
        calculator (ase.calculators.calculator.Calculator): The calculator to use for the energy calculation

    Returns:
        float: The interfacial energy of the interface.
    """
    substrate_slab = get_slab(interface, part="substrate") if substrate_slab is None else substrate_slab
    film_slab = get_slab(interface, part="film") if film_slab is None else film_slab

    build_configuration = interface.metadata["build"]["configuration"] if "build" in interface.metadata else {}
    try:
        substrate_bulk = (
            Material(build_configuration["substrate_configuration"]["bulk"])
            if substrate_bulk is None
            else substrate_bulk
        )
        film_bulk = Material(build_configuration["film_configuration"]["bulk"]) if film_bulk is None else film_bulk
    except KeyError:
        raise ValueError("The substrate and film bulk materials must be provided or defined in the interface metadata.")
    surface_energy_substrate = calculate_surface_energy(substrate_slab, substrate_bulk, calculator)
    surface_energy_film = calculate_surface_energy(film_slab, film_bulk, calculator)
    adhesion_energy = calculate_adhesion_energy(interface, substrate_slab, film_slab, calculator)
    return surface_energy_film + surface_energy_substrate - adhesion_energy
