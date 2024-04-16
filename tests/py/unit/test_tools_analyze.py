import numpy as np
from ase.build import bulk
from mat3ra.made.tools.analyze import calculate_average_interlayer_distance


def test_calculate_average_interlayer_distance():
    substrate = bulk("Si", cubic=True)
    film = bulk("Cu", cubic=True)
    interface_atoms = substrate + film
    interface_atoms.set_tags([1] * len(substrate) + [2] * len(film))
    distance = calculate_average_interlayer_distance(interface_atoms, 1, 2)
    assert np.isclose(distance, 4.0725)
