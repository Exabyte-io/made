import numpy as np
from ase.build import bulk, surface
from ase.calculators import emt
from mat3ra.made.tools.calculate import (
    calculate_total_energy,
    calculate_total_energy_per_atom,
    calculate_surface_energy,
    calculate_adhesion_energy,
)


def test_calculate_total_energy():
    atoms = bulk("C", cubic=True)
    calculator = emt.EMT()
    energy = calculate_total_energy(atoms, calculator)
    assert np.isclose(energy, 1.3612647524769237)


def test_calculate_total_energy_per_atom():
    atoms = bulk("C", cubic=True)
    calculator = emt.EMT()
    print(atoms.get_global_number_of_atoms())
    energy_per_atom = calculate_total_energy_per_atom(atoms, calculator)
    assert np.isclose(energy_per_atom, 0.1701580940596)


def test_calculate_surface_energy():
    atoms_slab = surface("C", (1, 1, 1), 3, vacuum=10)
    atoms_bulk = bulk("C", cubic=True)
    calculator = emt.EMT()
    surface_energy = calculate_surface_energy(atoms_slab, atoms_bulk, calculator)
    assert np.isclose(surface_energy, 0.148845)
