import numpy as np
from ase.build import bulk, surface
from mat3ra.made.tools.analyze import calculate_average_interlayer_distance
from mat3ra.made.tools.analyze import calculate_surface_area


def test_calculate_average_interlayer_distance():
    substrate = bulk("Si", cubic=True)
    film = bulk("Cu", cubic=True)
    interface_atoms = substrate + film
    interface_atoms.set_tags([1] * len(substrate) + [2] * len(film))
    distance = calculate_average_interlayer_distance(interface_atoms, 1, 2)
    assert np.isclose(distance, 4.0725)


def test_calculate_surface_area():
    atoms = bulk("Si", cubic=False)
    area = calculate_surface_area(atoms)
    assert np.isclose(area, 12.7673)
