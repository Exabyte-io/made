from types import SimpleNamespace
from typing import Final

import numpy as np
import pytest
from mat3ra.made.tools.analyze.interface import InterfaceAnalyzer
from mat3ra.made.tools.analyze.interface.commensurate import CommensurateLatticeInterfaceAnalyzer
from mat3ra.made.tools.build.slab.configurations import SlabConfiguration
from unit.fixtures.bulk import BULK_Ge_CONVENTIONAL, BULK_Si_CONVENTIONAL

from .fixtures.monolayer import GRAPHENE
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
    substrate_slab_config = SlabConfiguration.from_parameters(
        substrate.bulk_config, substrate.miller_indices, substrate.number_of_layers, vacuum=substrate.vacuum
    )
    film_slab_config = SlabConfiguration.from_parameters(
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

    # Test strain matrix calculation using vectors from substrate and film
    substrate_vectors = np.array(interface_analyzer.substrate_material.lattice.vector_arrays)
    film_vectors = np.array(interface_analyzer.film_material.lattice.vector_arrays)

    strain_matrix_3d = interface_analyzer.get_film_strain_matrix(substrate_vectors, film_vectors)

    matrix = [row.root for row in strain_matrix_3d.root]
    strain_matrix_2d = np.array(matrix)[:2, :2]

    # Verify that strain applied to film vectors yields substrate vectors
    film_2d = film_vectors[:2, :2]
    substrate_2d = substrate_vectors[:2, :2]
    film_strained = film_2d @ strain_matrix_2d

    assert np.allclose(film_strained, substrate_2d, atol=1e-10)


@pytest.mark.parametrize(
    "material_config, analyzer_params, expected_matches_len, expected_angle_range",
    [
        (
            GRAPHENE,
            {"target_angle": 13.0, "angle_tolerance": 0.5, "max_supercell_matrix_int": 5, "return_first_match": True},
            1,
            (12.5, 13.5),
        )
    ],
)
def test_commensurate_analyzer_functionality(
    material_config, analyzer_params, expected_matches_len, expected_angle_range
):
    slab_config = SlabConfiguration.from_parameters(
        material_config, miller_indices=(0, 0, 1), number_of_layers=1, vacuum=0.0
    )

    analyzer = CommensurateLatticeInterfaceAnalyzer(substrate_slab_configuration=slab_config, **analyzer_params)

    assert analyzer.substrate_material == analyzer.film_material == analyzer.material

    match_holders = analyzer.commensurate_lattice_match_holders
    assert len(match_holders) >= expected_matches_len

    for match in match_holders:
        assert expected_angle_range[0] <= match.angle <= expected_angle_range[1]
        assert expected_angle_range[0] <= match.actual_angle <= expected_angle_range[1]
        assert match.angle == match.actual_angle
        assert isinstance(match.xy_supercell_matrix_film, list)
        assert isinstance(match.xy_supercell_matrix_substrate, list)
        assert len(match.xy_supercell_matrix_film) == 2
        assert len(match.xy_supercell_matrix_substrate) == 2
        assert isinstance(match.match_id, int)
        assert match.match_id >= 0

    interface_configurations = analyzer.get_strained_configurations()
    assert len(interface_configurations) == len(match_holders)

    if len(match_holders) > 0:
        selected_config = analyzer.get_strained_configuration_by_match_id(0)
        assert selected_config.match_id == 0
        assert hasattr(selected_config, "substrate_configuration")
        assert hasattr(selected_config, "film_configuration")

        assert selected_config.substrate_configuration.stack_components == slab_config.stack_components
        assert selected_config.film_configuration.stack_components == slab_config.stack_components

        # Test invalid match ID
        with pytest.raises(ValueError, match="Match ID .* out of range"):
            analyzer.get_strained_configuration_by_match_id(999)

        # Test negative match ID
        with pytest.raises(ValueError, match="Match ID .* out of range"):
            analyzer.get_strained_configuration_by_match_id(-1)
