from types import SimpleNamespace
from typing import Final

import numpy as np
import pytest
from mat3ra.made.tools.analyze.interface import InterfaceAnalyzer, StrainedSlabConfigurationHolder, ZSLInterfaceAnalyzer
from mat3ra.made.tools.analyze.utils import get_film_strain_matrix
from mat3ra.made.tools.build.slab.helpers import create_slab_configuration
from unit.fixtures.bulk import BULK_Ge_CONVENTIONAL, BULK_Si_CONVENTIONAL

from .utils import assert_two_entities_deep_almost_equal

SUBSTRATE_SI_001: Final = SimpleNamespace(
    bulk_config=BULK_Si_CONVENTIONAL,
    miller_indices=(0, 0, 1),
    number_of_layers=2,
    vacuum=0.0,
)

FILM_GE_001: Final = SimpleNamespace(
    bulk_config=BULK_Ge_CONVENTIONAL,
    miller_indices=(0, 0, 1),
    number_of_layers=2,
    vacuum=0.0,
)

EXPECTED_PROPERTIES_SI_GE_001: Final = SimpleNamespace(
    substrate_strain_matrix=np.identity(3).tolist(),
    substrate_supercell_matrix=[[1.0, 0.0], [0.0, 1.0]],
    film_strain_matrix=[[0.9643, 0.0, 0.0], [0.0, 0.9643, 0.0], [0.0, 0.0, 1.0]],
    film_supercell_matrix=[[1.0, 0.0], [0.0, 1.0]],
)

TEST_CASES = [(SUBSTRATE_SI_001, FILM_GE_001, EXPECTED_PROPERTIES_SI_GE_001)]


@pytest.mark.parametrize("substrate, film, expected", TEST_CASES)
def test_interface_analyzer(substrate, film, expected):
    substrate_slab_config = create_slab_configuration(
        substrate.bulk_config, substrate.miller_indices, substrate.number_of_layers, vacuum=substrate.vacuum
    )
    film_slab_config = create_slab_configuration(
        film.bulk_config, film.miller_indices, film.number_of_layers, vacuum=film.vacuum
    )

    interface_analyzer = InterfaceAnalyzer(
        substrate_slab_configuration=substrate_slab_config,
        film_slab_configuration=film_slab_config,
    )

    assert_two_entities_deep_almost_equal(
        interface_analyzer.get_substrate_strain_matrix().root, expected.substrate_strain_matrix
    )
    assert_two_entities_deep_almost_equal(
        interface_analyzer.film_strained_configuration.strain_matrix.root, expected.film_strain_matrix, atol=1e-4
    )


@pytest.mark.parametrize(
    "substrate_vectors, film_vectors, expected_strain_2d",
    [
        (
            np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
            np.array([[0.9, 0.0, 0.0], [0.0, 1.2, 0.0], [0.0, 0.0, 1.0]]),
            np.array([[1.0 / 0.9, 0.0], [0.0, 1.0 / 1.2]]),
        ),
        (
            np.array([[2.0, 0.0, 0.0], [0.0, 6.0, 0.0], [0.0, 0.0, 3.0]]),
            np.array([[1.8, 0.0, 0.0], [3.6, 7.2, 0.0], [0.0, 0.0, 2.5]]),
            np.array([[2.0 / 1.8, 0.0], [-2.0 / 3.6, 6.0 / 7.2]]),
        ),
    ],
)
def test_get_film_strain_matrix(substrate_vectors, film_vectors, expected_strain_2d):
    strain_matrix_3d = get_film_strain_matrix(substrate_vectors, film_vectors)
    strain_matrix_2d = strain_matrix_3d[:2, :2]

    assert np.allclose(strain_matrix_2d, expected_strain_2d, atol=1e-10)

    # Verify that strain applied to film vectors yields substrate vectors
    film_2d = film_vectors[:2, :2]
    substrate_2d = substrate_vectors[:2, :2]
    film_strained = film_2d @ strain_matrix_2d

    assert np.allclose(film_strained, substrate_2d, atol=1e-10)
