import numpy as np
from ase.build import bulk
from ase.calculators import emt
from mat3ra.made.tools.calculate import calculate_total_energy


def test_calculate_total_energy():
    atoms = bulk("C", cubic=True)
    calculator = emt.EMT()
    energy = calculate_total_energy(atoms, calculator)
    assert np.isclose(energy, 1.3612647524769237)
