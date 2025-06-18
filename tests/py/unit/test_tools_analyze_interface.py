import numpy as np
import pytest
from mat3ra.utils.assertion import assert_deep_almost_equal

from mat3ra.made.tools.analyze.interface import InterfaceAnalyzer
from unit.fixtures.bulk import BULK_Ge_CONVENTIONAL, BULK_Si_CONVENTIONAL
from .helpers import get_slab_configuration, calculate_expected_film_strain_matrix


@pytest.fixture(name="interface_analyzer")
def get_interface_analyzer() -> InterfaceAnalyzer:
    substrate_slab_config = get_slab_configuration(BULK_Si_CONVENTIONAL, (0, 0, 1), 2, vacuum=0.0)
    film_slab_config = get_slab_configuration(BULK_Ge_CONVENTIONAL, (0, 0, 1), 2, vacuum=0.0)
    return InterfaceAnalyzer(
        substrate_slab_configuration=substrate_slab_config,
        film_slab_configuration=film_slab_config,
    )


def test_substrate_strain_matrix(interface_analyzer: InterfaceAnalyzer):
    assert_deep_almost_equal(interface_analyzer.substrate_strain_matrix.root, np.identity(3).tolist())


def test_substrate_supercell_matrix(interface_analyzer: InterfaceAnalyzer):
    assert_deep_almost_equal(interface_analyzer.substrate_supercell_matrix.root, [[1.0, 0.0], [0.0, 1.0]])


def test_film_supercell_matrix(interface_analyzer: InterfaceAnalyzer):
    assert_deep_almost_equal(interface_analyzer.film_supercell_matrix.root, [[1.0, 0.0], [0.0, 1.0]])


def test_film_strain_matrix(interface_analyzer: InterfaceAnalyzer):
    film_material = interface_analyzer.film_material
    substrate_material = interface_analyzer.substrate_material
    expected_film_strain_matrix = calculate_expected_film_strain_matrix(film_material, substrate_material)
    assert_deep_almost_equal(
        interface_analyzer.film_strain_matrix.root,
        expected_film_strain_matrix.tolist(),
    )
