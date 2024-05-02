import numpy as np
from ase.build import bulk
from mat3ra.made.tools.analyze import get_average_interlayer_distance, get_surface_area

from .fixtures import INTERFACE_ATOMS


def test_calculate_average_interlayer_distance():
    distance = get_average_interlayer_distance(INTERFACE_ATOMS, 1, 2)
    assert np.isclose(distance, 4.0725)


def test_calculate_surface_area():
    atoms = bulk("Si", cubic=False)
    area = get_surface_area(atoms)
    assert np.isclose(area, 12.7673)
