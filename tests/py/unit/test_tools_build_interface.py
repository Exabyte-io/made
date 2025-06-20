import platform
from types import SimpleNamespace
from typing import Final

import pytest
from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.interface import InterfaceAnalyzer, ZSLInterfaceAnalyzer
from mat3ra.made.tools.build.interface import (
    InterfaceBuilder,
    InterfaceConfiguration,
    create_interfaces,
)
from mat3ra.made.tools.build.interface.builders import (
    CommensurateLatticeTwistedInterfaceBuilder,
    CommensurateLatticeTwistedInterfaceBuilderParameters,
    NanoRibbonTwistedInterfaceBuilder,
    NanoRibbonTwistedInterfaceConfiguration,
    TwistedInterfaceConfiguration,
)
from mat3ra.made.tools.build.slab.helpers import create_slab_configuration
from mat3ra.utils import assertion as assertion_utils
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


def test_zsl_interface_builder():
    """Test creating Si/Ge interface using ZSL approach."""
    # Create slab configurations for Si (substrate) and Ge (film)
    substrate_slab_config = create_slab_configuration(BULK_Si_CONVENTIONAL, (0, 0, 1), 3, vacuum=10.0)
    film_slab_config = create_slab_configuration(BULK_Ge_CONVENTIONAL, (0, 0, 1), 3, vacuum=10.0)

    # Use ZSLInterfaceAnalyzer to get strained slab configurations
    analyzer = ZSLInterfaceAnalyzer(
        substrate_slab_configuration=substrate_slab_config,
        film_slab_configuration=film_slab_config,
        max_area=50.0,  # Keep area reasonable for testing
        max_area_ratio_tol=0.1,
        max_length_tol=0.05,
        max_angle_tol=0.02,
    )

    # Get strained configurations with metadata
    configs_with_metadata = analyzer.get_strained_slab_configurations_with_metadata()

    # Should find at least one ZSL match
    assert len(configs_with_metadata) > 0, "No ZSL matches found"

    # Select the configuration with lowest strain
    selected_config = configs_with_metadata[0]

    # Verify strain information is present
    assert "substrate_strain" in selected_config.strain_info
    assert "film_strain" in selected_config.strain_info
    assert selected_config.strain_info["substrate_strain"] >= 0
    assert selected_config.strain_info["film_strain"] >= 0

    # Create interface configuration using the strained slab configs
    interface_config = InterfaceConfiguration(
        stack_components=[
            selected_config.substrate_config,
            selected_config.film_config,
        ]
    )

    # Use regular InterfaceBuilder to create the interface
    builder = InterfaceBuilder()
    interface = builder.get_material(interface_config)

    # Should contain both Si and Ge
    elements = set(interface.basis.elements.values)
    assert "Si" in elements, "Interface should contain Si atoms"
    assert "Ge" in elements, "Interface should contain Ge atoms"


@pytest.mark.skip(reason="Fixtures are commented out. To be fixed in epic-7623")
def test_create_interfaces():
    # interfaces = create_interfaces(matched_interfaces_builder, interface_configuration)

    # assert len(interfaces) == EXPECTED_NUMBER_OF_INTERFACES
    # assert interfaces[0].name is not None
    pass


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
    "material_config, config_params, builder_params_dict,"
    + " expected_interfaces_len, expected_cell_vectors, expected_angle",
    [
        (
            GRAPHENE,
            {"twist_angle": 13, "distance_z": 3.0},
            {"max_supercell_matrix_int": 5, "angle_tolerance": 0.5, "return_first_match": True},
            1,
            [[10.754672133, 0.0, 0.0], [5.377336066500001, 9.313819276550575, 0.0], [0.0, 0.0, 20.0]],
            13.174,
        ),
    ],
)
def test_create_commensurate_supercell_twisted_interface(
    material_config,
    config_params,
    builder_params_dict,
    expected_interfaces_len,
    expected_cell_vectors,
    expected_angle,
):
    film = Material.create(material_config)
    substrate = Material.create(material_config)
    config = TwistedInterfaceConfiguration(film=film, substrate=substrate, **config_params)
    params = CommensurateLatticeTwistedInterfaceBuilderParameters(**builder_params_dict)
    builder = CommensurateLatticeTwistedInterfaceBuilder(build_parameters=params)
    interfaces = builder.get_materials(config, post_process_parameters=config)
    assert len(interfaces) == expected_interfaces_len
    interface = interfaces[0]
    assertion_utils.assert_deep_almost_equal(expected_cell_vectors, interface.basis.cell.vector_arrays)
    assert interface.metadata["build"]["configuration"]["actual_twist_angle"] == expected_angle
