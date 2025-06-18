import numpy as np
import pytest
from mat3ra.made.tools.analyze.interface import InterfaceAnalyzer
from mat3ra.utils.assertion import assert_deep_almost_equal
from unit.fixtures.bulk import BULK_Ge_CONVENTIONAL, BULK_Si_CONVENTIONAL

from .helpers import get_slab_configuration

TEST_CASES = [
    (
        BULK_Si_CONVENTIONAL,
        (0, 0, 1),
        2,
        0.0,
        BULK_Ge_CONVENTIONAL,
        (0, 0, 1),
        2,
        0.0,
        np.identity(3).tolist(),
        [[1.0, 0.0], [0.0, 1.0]],
        [[0.96, 0.0, 0.0], [0.0, 0.96, 0.0], [0.0, 0.0, 1.0]],
        [[1.0, 0.0], [0.0, 1.0]],
    )
]


@pytest.mark.parametrize(
    "substrate_bulk_config, substrate_miller_indices, substrate_number_of_layers, substrate_vacuum, "
    + " film_bulk_config, film_miller_indices, film_number_of_layers, film_vacuum, "
    + " expected_substrate_strain_matrix, expected_substrate_supercell_matrix,"
    + " expected_film_strain_matrix, expected_film_supercell_matrix",
    TEST_CASES,
)
def test_interface_analyzer(
    substrate_bulk_config,
    substrate_miller_indices,
    substrate_number_of_layers,
    substrate_vacuum,
    film_bulk_config,
    film_miller_indices,
    film_number_of_layers,
    film_vacuum,
    expected_substrate_strain_matrix,
    expected_substrate_supercell_matrix,
    expected_film_strain_matrix,
    expected_film_supercell_matrix,
):
    substrate_slab_config = get_slab_configuration(
        substrate_bulk_config, substrate_miller_indices, substrate_number_of_layers, vacuum=substrate_vacuum
    )
    film_slab_config = get_slab_configuration(
        film_bulk_config, film_miller_indices, film_number_of_layers, vacuum=film_vacuum
    )
    interface_analyzer = InterfaceAnalyzer(
        substrate_slab_configuration=substrate_slab_config,
        film_slab_configuration=film_slab_config,
    )

    assert_deep_almost_equal(interface_analyzer.substrate_strain_matrix.root, expected_substrate_strain_matrix)
    assert_deep_almost_equal(
        interface_analyzer.substrate_supercell_matrix.root,
        expected_substrate_supercell_matrix,
    )

    assert_deep_almost_equal(
        interface_analyzer.film_strain_matrix.root,
        expected_film_strain_matrix,
    )
    assert_deep_almost_equal(interface_analyzer.film_supercell_matrix.root, expected_film_supercell_matrix)
