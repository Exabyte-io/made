import platform
from types import SimpleNamespace

import pytest
from mat3ra.code.array_with_ids import ArrayWithIds
from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum

from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.interface.commensurate import CommensurateLatticeInterfaceAnalyzer
from mat3ra.made.tools.analyze.interface.simple import InterfaceAnalyzer
from mat3ra.made.tools.analyze.interface.zsl import ZSLInterfaceAnalyzer
from mat3ra.made.tools.build.interface.builders import InterfaceBuilder, InterfaceConfiguration
from mat3ra.made.tools.build.interface.helpers import (
    create_simple_interface_between_slabs,
    create_twisted_interface,
    create_zsl_interface,
)
from mat3ra.made.tools.build.nanoribbon import create_nanoribbon
from mat3ra.made.tools.build.slab.builders import SlabBuilder
from mat3ra.made.tools.build.slab.configurations import SlabConfiguration
from unit.fixtures.bulk import BULK_Ge_CONVENTIONAL, BULK_Si_CONVENTIONAL
from .fixtures.interface.commensurate import INTERFACE_GRAPHENE_GRAPHENE_X, INTERFACE_GRAPHENE_GRAPHENE_Z
from .fixtures.interface.simple import INTERFACE_Si_001_Ge_001  # type: ignore
from .fixtures.interface.twisted_nanoribbons import TWISTED_INTERFACE_GRAPHENE_GRAPHENE_60
from .fixtures.monolayer import GRAPHENE
from .utils import assert_two_entities_deep_almost_equal

SIMPLE_INTERFACE_BUILDER_TEST_CASES = [
    (
        SimpleNamespace(
            bulk_config=BULK_Si_CONVENTIONAL,
            miller_indices=(0, 0, 1),
            number_of_layers=2,
            vacuum=0.0,
        ),
        SimpleNamespace(
            bulk_config=BULK_Ge_CONVENTIONAL,
            miller_indices=(0, 0, 1),
            number_of_layers=2,
            vacuum=0.0,
        ),
        INTERFACE_Si_001_Ge_001,
    )
]


MAX_AREA = 100
# pymatgen `2023.6.23` supporting py3.8 returns 1 interface instead of 2
EXPECTED_NUMBER_OF_INTERFACES = 1 if platform.python_version().startswith("3.8") else 1
interface_configuration = None


@pytest.mark.parametrize("substrate, film, expected_interface", SIMPLE_INTERFACE_BUILDER_TEST_CASES)
def test_simple_interface_builder(substrate, film, expected_interface):
    builder = InterfaceBuilder()
    substrate_slab_config = SlabConfiguration.from_parameters(
        substrate.bulk_config, substrate.miller_indices, substrate.number_of_layers, vacuum=substrate.vacuum
    )
    film_slab_config = SlabConfiguration.from_parameters(
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
    substrate_slab_config = SlabConfiguration.from_parameters(
        substrate.bulk_config, substrate.miller_indices, substrate.number_of_layers, vacuum=substrate.vacuum
    )
    film_slab_config = SlabConfiguration.from_parameters(
        film.bulk_config, film.miller_indices, film.number_of_layers, vacuum=film.vacuum
    )

    # Use ZSLInterfaceAnalyzer to get strained slab configurations
    analyzer = ZSLInterfaceAnalyzer(
        substrate_slab_configuration=substrate_slab_config,
        film_slab_configuration=film_slab_config,
        max_area=MAX_AREA,
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
    interface.metadata.build = []
    expected_interface["metadata"].pop("build", None)

    assert_two_entities_deep_almost_equal(interface, expected_interface)


@pytest.mark.parametrize("substrate, film, expected_interface", SIMPLE_INTERFACE_BUILDER_TEST_CASES)
def test_create_simple_interface_between_slabs(substrate, film, expected_interface):
    substrate_slab_config = SlabConfiguration.from_parameters(
        material_or_dict=substrate.bulk_config,
        miller_indices=substrate.miller_indices,
        number_of_layers=substrate.number_of_layers,
        vacuum=0,
        termination_formula=None,
    )
    film_slab_config = SlabConfiguration.from_parameters(
        material_or_dict=film.bulk_config,
        miller_indices=film.miller_indices,
        number_of_layers=film.number_of_layers,
        vacuum=0,
        termination_formula=None,
    )

    substrate_slab = SlabBuilder().get_material(substrate_slab_config)
    film_slab = SlabBuilder().get_material(film_slab_config)

    interface = create_simple_interface_between_slabs(
        substrate_slab=substrate_slab,
        film_slab=film_slab,
        gap=None,
        vacuum=0.0,
        xy_shift=[0, 0],
    )
    interface.metadata.build = []
    expected_interface["metadata"].pop("build", None)
    assert_two_entities_deep_almost_equal(interface, expected_interface)


@pytest.mark.parametrize("substrate, film, expected_interface", SIMPLE_INTERFACE_BUILDER_TEST_CASES)
def test_create_zsl_interface(substrate, film, expected_interface):
    interface = create_zsl_interface(
        substrate_crystal=substrate.bulk_config,
        film_crystal=film.bulk_config,
        substrate_miller_indices=substrate.miller_indices,
        film_miller_indices=film.miller_indices,
        substrate_number_of_layers=substrate.number_of_layers,
        film_number_of_layers=film.number_of_layers,
        substrate_termination_formula=None,
        film_termination_formula=None,
        gap=None,
        vacuum=0.0,
        xy_shift=[0, 0],
        max_area=MAX_AREA,
        max_area_ratio_tol=0.1,
        max_length_tol=0.05,
        max_angle_tol=0.02,
    )
    interface.metadata.build = []
    expected_interface["metadata"].pop("build", None)
    assert_two_entities_deep_almost_equal(interface, expected_interface)


@pytest.mark.parametrize(
    "material_config, config_params, expected_interface",
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
            TWISTED_INTERFACE_GRAPHENE_GRAPHENE_60,
        ),
    ],
)
def test_create_twisted_interface(material_config, config_params, expected_interface):
    material = Material.create(material_config)

    nanoribbon1 = create_nanoribbon(
        material=material,
        width=3,
        length=3,
        vacuum_width=15.0,
        vacuum_length=15.0,
    )

    nanoribbon2 = create_nanoribbon(
        material=material,
        width=3,
        length=3,
        vacuum_width=15.0,
        vacuum_length=15.0,
    )

    interface = create_twisted_interface(
        material1=nanoribbon1,
        material2=nanoribbon2,
        angle=config_params["twist_angle"],
        vacuum_x=config_params["vacuum_x"],
        vacuum_y=config_params["vacuum_y"],
        gap=config_params["distance_z"],
    )

    assert_two_entities_deep_almost_equal(interface, expected_interface)


@pytest.mark.parametrize(
    "material_config, analyzer_params, direction, gaps, expected_interface",
    [
        (
            GRAPHENE,
            {"target_angle": 13.0, "angle_tolerance": 0.5, "max_supercell_matrix_int": 5, "return_first_match": True},
            AxisEnum.z,
            [3.0, 3.0],
            INTERFACE_GRAPHENE_GRAPHENE_Z,
        ),
        (
            GRAPHENE,
            {"target_angle": 13.0, "angle_tolerance": 0.5, "max_supercell_matrix_int": 5, "return_first_match": True},
            AxisEnum.x,
            [],  # testing no gaps
            INTERFACE_GRAPHENE_GRAPHENE_X,
        ),
    ],
)
def test_commensurate_interface_creation(material_config, analyzer_params, direction, gaps, expected_interface):
    slab_config = SlabConfiguration.from_parameters(
        material_config, miller_indices=(0, 0, 1), number_of_layers=1, vacuum=0.0
    )

    analyzer = CommensurateLatticeInterfaceAnalyzer(substrate_slab_configuration=slab_config, **analyzer_params)

    match_holders = analyzer.commensurate_lattice_match_holders

    if len(match_holders) > 0:
        selected_config = analyzer.get_strained_configuration_by_match_id(0)
        interface_config = InterfaceConfiguration(
            stack_components=[selected_config.substrate_configuration, selected_config.film_configuration],
            gaps=ArrayWithIds.from_values(values=gaps),
            direction=direction,
        )

        builder = InterfaceBuilder()
        interface = builder.get_material(interface_config)
        interface.metadata.build = []

        assert_two_entities_deep_almost_equal(interface, expected_interface, atol=1e-5)
