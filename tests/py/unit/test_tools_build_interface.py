import platform
from types import SimpleNamespace
from typing import Final

import pytest
from mat3ra.utils import assertion as assertion_utils

from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.interface import InterfaceAnalyzer, ZSLInterfaceAnalyzer
from mat3ra.made.tools.analyze.interface.commensurate import CommensurateLatticeTwistedInterfaceAnalyzer
from mat3ra.made.tools.build.interface import InterfaceBuilder, InterfaceConfiguration, create_interface
from mat3ra.made.tools.build.interface.builders import (
    NanoRibbonTwistedInterfaceBuilder,
    NanoRibbonTwistedInterfaceConfiguration,
)
from mat3ra.made.tools.build.slab.helpers import create_slab_configuration
from unit.fixtures.bulk import BULK_Ge_CONVENTIONAL, BULK_Si_CONVENTIONAL
from .fixtures.interface.simple import INTERFACE_Si_001_Ge_001  # type: ignore
from .fixtures.monolayer import GRAPHENE
from .utils import assert_two_entities_deep_almost_equal

PARAMETERS_SLAB_Si_001: Final = SimpleNamespace(
    bulk_config=BULK_Si_CONVENTIONAL,
    miller_indices=(0, 0, 1),
    number_of_layers=2,
    vacuum=0.0,
)

PARAMETERS_SLAB_Ge_001: Final = SimpleNamespace(
    bulk_config=BULK_Ge_CONVENTIONAL,
    miller_indices=(0, 0, 1),
    number_of_layers=2,
    vacuum=0.0,
)


SIMPLE_INTERFACE_BUILDER_TEST_CASES = [(PARAMETERS_SLAB_Si_001, PARAMETERS_SLAB_Ge_001, INTERFACE_Si_001_Ge_001)]


MAX_AREA = 100
# pymatgen `2023.6.23` supporting py3.8 returns 1 interface instead of 2
EXPECTED_NUMBER_OF_INTERFACES = 1 if platform.python_version().startswith("3.8") else 1
interface_configuration = None


@pytest.mark.parametrize("substrate, film, expected_interface", SIMPLE_INTERFACE_BUILDER_TEST_CASES)
def test_simple_interface_builder(substrate, film, expected_interface):
    builder = InterfaceBuilder()
    substrate_slab_config = create_slab_configuration(
        substrate.bulk_config, substrate.miller_indices, substrate.number_of_layers, vacuum=substrate.vacuum
    )
    film_slab_config = create_slab_configuration(
        film.bulk_config, film.miller_indices, film.number_of_layers, vacuum=film.vacuum
    )

    analyzer = InterfaceAnalyzer(
        substrate_slab_configuration=substrate_slab_config,
        film_slab_configuration=film_slab_config,
    )

    film_configuration = analyzer.film_strained_configuration
    substrate_configuration = analyzer.substrate_strained_configuration
    vacuum_configuration = film_slab_config.vacuum_configuration

    config = InterfaceConfiguration(
        stack_components=[substrate_configuration, film_configuration, vacuum_configuration],
    )
    interface = builder.get_material(config)
    assert_two_entities_deep_almost_equal(interface, expected_interface)


@pytest.mark.parametrize("substrate, film, expected_interface", SIMPLE_INTERFACE_BUILDER_TEST_CASES)
def test_zsl_interface_builder(substrate, film, expected_interface):
    """Test creating Si/Ge interface using ZSL approach."""
    substrate_slab_config = create_slab_configuration(
        substrate.bulk_config, substrate.miller_indices, substrate.number_of_layers, vacuum=substrate.vacuum
    )
    film_slab_config = create_slab_configuration(
        film.bulk_config, film.miller_indices, film.number_of_layers, vacuum=film.vacuum
    )

    # Use ZSLInterfaceAnalyzer to get strained slab configurations
    analyzer = ZSLInterfaceAnalyzer(
        substrate_slab_configuration=substrate_slab_config,
        film_slab_configuration=film_slab_config,
        max_area=50.0,
        max_area_ratio_tol=0.1,
        max_length_tol=0.05,
        max_angle_tol=0.02,
    )

    interface_configurations = analyzer.get_strained_configurations()

    selected_config = interface_configurations[0]
    interface_config = InterfaceConfiguration(
        stack_components=[selected_config.substrate_configuration, selected_config.film_configuration]
    )

    builder = InterfaceBuilder()
    interface = builder.get_material(interface_config)

    # remove metadata
    interface.metadata.pop("build", None)
    expected_interface["metadata"].pop("build", None)

    assert_two_entities_deep_almost_equal(interface, expected_interface)


@pytest.mark.parametrize("substrate, film, expected_interface", SIMPLE_INTERFACE_BUILDER_TEST_CASES)
def test_create_interface(substrate, film, expected_interface):
    substrate_slab_config = create_slab_configuration(
        substrate.bulk_config, substrate.miller_indices, substrate.number_of_layers, vacuum=substrate.vacuum
    )
    film_slab_config = create_slab_configuration(
        film.bulk_config, film.miller_indices, film.number_of_layers, vacuum=film.vacuum
    )

    analyzer = InterfaceAnalyzer(
        substrate_slab_configuration=substrate_slab_config,
        film_slab_configuration=film_slab_config,
    )

    film_configuration = analyzer.film_strained_configuration
    substrate_configuration = analyzer.substrate_strained_configuration
    vacuum_configuration = film_slab_config.vacuum_configuration

    configuration = InterfaceConfiguration(
        stack_components=[substrate_configuration, film_configuration, vacuum_configuration],
    )

    interface = create_interface(configuration)
    interface.metadata.pop("build", None)
    expected_interface["metadata"].pop("build", None)


@pytest.mark.parametrize(
    "material_config, config_params, expected_cell_vectors, expected_coordinate_checks",
    [
        (
            GRAPHENE,
            {
                "twist_angle": 60,
                "distance_z": 3.0,
                "ribbon_width": 3,
                "ribbon_length": 5,
                "vacuum_x": 2.0,
                "vacuum_y": 2.0,
            },
            [[15.102811, 0.0, 0.0], [0.0, 16.108175208, 0.0], [0.0, 0.0, 20.0]],
            {42: [0.704207885, 0.522108183, 0.65]},
        ),
    ],
)
def test_create_twisted_nanoribbon_interface(
    material_config, config_params, expected_cell_vectors, expected_coordinate_checks
):
    film = Material.create(material_config)
    configuration = NanoRibbonTwistedInterfaceConfiguration(film=film, substrate=film, **config_params)

    builder = NanoRibbonTwistedInterfaceBuilder()
    interface = builder.get_material(configuration)

    assertion_utils.assert_deep_almost_equal(expected_cell_vectors, interface.lattice.vector_arrays)
    for index, expected_coordinate in expected_coordinate_checks.items():
        assertion_utils.assert_deep_almost_equal(expected_coordinate, interface.basis.coordinates.values[index])


@pytest.mark.parametrize(
    "material_config, analyzer_params, expected_matches_len, expected_angle_range",
    [
        (
            GRAPHENE,
            {"target_angle": 13.0, "angle_tolerance": 0.5, "max_supercell_matrix_int": 5, "return_first_match": True},
            1,
            (12.5, 13.5),
        ),
    ],
)
def test_commensurate_lattice_twisted_interface_analyzer(
    material_config, analyzer_params, expected_matches_len, expected_angle_range
):

    # Create slab configuration (both film and substrate use the same material for twisted bilayers)
    slab_config = create_slab_configuration(material_config, miller_indices=(0, 0, 1), number_of_layers=1, vacuum=0.0)

    analyzer = CommensurateLatticeTwistedInterfaceAnalyzer(substrate_slab_configuration=slab_config, **analyzer_params)

    match_holders = analyzer.commensurate_match_holders
    assert len(match_holders) >= expected_matches_len

    for match in match_holders:
        assert expected_angle_range[0] <= match.angle <= expected_angle_range[1]
        assert isinstance(match.xy_supercell_matrix_film, list)
        assert isinstance(match.xy_supercell_matrix_substrate, list)
        assert len(match.xy_supercell_matrix_film) == 2
        assert len(match.xy_supercell_matrix_substrate) == 2

    interface_configurations = analyzer.get_strained_configurations()
    assert len(interface_configurations) == len(match_holders)

    if len(match_holders) > 0:
        selected_config = analyzer.get_strained_configuration_by_match_id(0)

        # Test that supercell matrices are properly set
        assert (
            selected_config.substrate_configuration.xy_supercell_matrix
            == match_holders[0].xy_supercell_matrix_substrate
        )
        assert selected_config.film_configuration.xy_supercell_matrix == match_holders[0].xy_supercell_matrix_film

        # Test building interface from analyzer configurations
        interface_config = InterfaceConfiguration(
            stack_components=[selected_config.substrate_configuration, selected_config.film_configuration]
        )

        builder = InterfaceBuilder()
        interface = builder.get_material(interface_config)

        # Check against expected interface
        assert isinstance(interface, Material)
