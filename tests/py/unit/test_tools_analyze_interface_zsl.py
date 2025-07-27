from types import SimpleNamespace
from typing import Final

import numpy as np
import pytest
from mat3ra.standata.materials import Materials
from mat3ra.utils.matrix import convert_2x2_to_3x3

from mat3ra.made.tools.analyze.interface.zsl import ZSLInterfaceAnalyzer
from mat3ra.made.tools.build.slab.builders import SlabStrainedSupercellBuilder
from mat3ra.made.tools.build.slab.configurations import SlabConfiguration
from unit.fixtures.bulk import BULK_Ge_CONVENTIONAL, BULK_Si_CONVENTIONAL
from .fixtures.monolayer import GRAPHENE

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


@pytest.mark.parametrize(
    "substrate, film, zsl_params, expected_matches_min",
    [
        (
            SUBSTRATE_SI_001,
            FILM_GE_001,
            {"max_area": 350.0, "max_area_ratio_tol": 0.09, "max_length_tol": 0.03, "max_angle_tol": 0.01},
            1,
        ),
        (
            SimpleNamespace(
                bulk_config=Materials.get_by_name_first_match("Nickel"),
                miller_indices=(0, 0, 1),
                number_of_layers=2,
                vacuum=0.0,
            ),
            SimpleNamespace(
                bulk_config=GRAPHENE,
                miller_indices=(0, 0, 1),
                number_of_layers=1,
                vacuum=0.0,
            ),
            {"max_area": 120.0, "max_area_ratio_tol": 0.09, "max_length_tol": 0.03, "max_angle_tol": 0.01},
            1,
        ),
    ],
)
def test_zsl_interface_analyzer(substrate, film, zsl_params, expected_matches_min):
    substrate_slab_config = SlabConfiguration.from_parameters(
        substrate.bulk_config, substrate.miller_indices, substrate.number_of_layers, vacuum=0.0
    )
    film_slab_config = SlabConfiguration.from_parameters(
        film.bulk_config, film.miller_indices, film.number_of_layers, vacuum=0.0
    )

    analyzer = ZSLInterfaceAnalyzer(
        substrate_slab_configuration=substrate_slab_config,
        film_slab_configuration=film_slab_config,
        **zsl_params,
    )

    # Test ZSL match generation
    match_holders = analyzer.zsl_match_holders
    assert len(match_holders) >= expected_matches_min

    match = analyzer.get_strained_configuration_by_match_id(0)

    sub_config = match.substrate_configuration
    film_config = match.film_configuration

    film_vectors = np.array(analyzer.film_slab.lattice.vector_arrays)
    substrate_vectors = np.array(analyzer.substrate_slab.lattice.vector_arrays)

    film_sl_vectors = (
        np.array(convert_2x2_to_3x3([row.root for row in film_config.xy_supercell_matrix]))
        @ film_vectors
        @ np.array([row.root for row in film_config.strain_matrix.root])
    )
    substrate_sl_vectors = (
        np.array(convert_2x2_to_3x3([row.root for row in sub_config.xy_supercell_matrix])) @ substrate_vectors
    )

    substrate_material = SlabStrainedSupercellBuilder().get_material(sub_config)
    film_material = SlabStrainedSupercellBuilder().get_material(film_config)

    assert np.allclose(film_sl_vectors[0:2], substrate_sl_vectors[0:2], atol=1e-4)

    assert np.allclose(substrate_material.lattice.vector_arrays[0:2], substrate_sl_vectors[0:2], atol=1e-4)
    assert np.allclose(film_material.lattice.vector_arrays[0:2], film_sl_vectors[0:2], atol=1e-4)

    assert np.isclose(substrate_material.lattice.a, film_material.lattice.a, atol=1e-4)
    assert np.isclose(substrate_material.lattice.b, film_material.lattice.b, atol=1e-4)


def test_zsl_interface_analyzer_sort_by_strain_then_area():
    substrate_slab_config = SlabConfiguration.from_parameters(
        BULK_Si_CONVENTIONAL, (0, 0, 1), number_of_layers=2, vacuum=0.0
    )
    film_slab_config = SlabConfiguration.from_parameters(
        BULK_Ge_CONVENTIONAL, (0, 0, 1), number_of_layers=2, vacuum=0.0
    )

    analyzer = ZSLInterfaceAnalyzer(
        substrate_slab_configuration=substrate_slab_config,
        film_slab_configuration=film_slab_config,
        max_area=350.0,
        max_area_ratio_tol=0.09,
        max_length_tol=0.03,
        max_angle_tol=0.01,
    )

    match_holders = analyzer.zsl_match_holders

    # Sort matches by strain percentage and then by area
    sorted_matches = sorted(match_holders, key=lambda x: (x.total_strain_percentage, x.match_area))

    assert len(sorted_matches) == len(match_holders)
    assert sorted_matches[0].total_strain_percentage <= sorted_matches[-1].total_strain_percentage
