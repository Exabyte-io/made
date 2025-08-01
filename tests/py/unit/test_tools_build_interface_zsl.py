from types import SimpleNamespace

import pytest
from mat3ra.made.tools.analyze.interface.zsl import ZSLInterfaceAnalyzer
from mat3ra.made.tools.build.interface.builders import InterfaceBuilder, InterfaceConfiguration
from mat3ra.made.tools.build.interface.helpers import create_zsl_interface, create_zsl_interface_between_slabs
from mat3ra.made.tools.build.slab.builders import SlabBuilder
from mat3ra.made.tools.build.slab.configurations import SlabConfiguration
from mat3ra.made.tools.build.vacuum.configuration import VacuumConfiguration
from mat3ra.standata.materials import Materials

from .fixtures.interface.gr_ni_111_top_hcp import (
    GRAPHENE_NICKEL_INTERFACE_TOP_HCP,
    GRAPHENE_NICKEL_INTERFACE_TOP_HCP_GH_WF,
)
from .fixtures.monolayer import GRAPHENE
from .utils import TestPlatform, assert_two_entities_deep_almost_equal, get_platform_specific_value

GRAPHENE_NICKEL_TEST_CASE = (
    SimpleNamespace(
        bulk_config=Materials.get_by_name_first_match("Nickel"),
        miller_indices=(1, 1, 1),
        number_of_layers=3,
        vacuum=0.0,
    ),
    SimpleNamespace(
        bulk_config=GRAPHENE,
        miller_indices=(0, 0, 1),
        number_of_layers=1,
        vacuum=0.0,
    ),
    3.0,  # gap between graphene and nickel
    10.0,  # vacuum
    90,  # max area
    {
        TestPlatform.DARWIN: GRAPHENE_NICKEL_INTERFACE_TOP_HCP,
        TestPlatform.OTHER: GRAPHENE_NICKEL_INTERFACE_TOP_HCP_GH_WF,
    },
)

MAX_AREA_RATIO_TOL = 0.09
MAX_LENGTH_TOL = 0.05
MAX_ANGLE_TOL = 0.02


@pytest.mark.parametrize(
    "substrate, film,gap, vacuum, max_area, expected_interface",
    [GRAPHENE_NICKEL_TEST_CASE],
)
def test_zsl_interface_builder(substrate, film, gap, vacuum, max_area, expected_interface):
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
        max_area=max_area,
        max_area_ratio_tol=MAX_AREA_RATIO_TOL,
        max_length_tol=MAX_LENGTH_TOL,
        max_angle_tol=MAX_ANGLE_TOL,
        reduce_result_cell=False,
    )

    interface_configurations = analyzer.get_strained_configurations()

    selected_config = interface_configurations[0]
    vacuum_configuration = VacuumConfiguration(size=vacuum)
    interface_config = InterfaceConfiguration(
        stack_components=[
            selected_config.substrate_configuration,
            selected_config.film_configuration,
            vacuum_configuration,
        ],
        gaps=[gap, gap],
    )

    builder = InterfaceBuilder()
    interface = builder.get_material(interface_config)

    # remove metadata
    interface.metadata.build = []
    expected_interface = get_platform_specific_value(expected_interface)
    print(interface.to_dict())
    assert_two_entities_deep_almost_equal(interface, expected_interface)


@pytest.mark.parametrize("substrate, film,gap, vacuum, max_area,  expected_interface", [GRAPHENE_NICKEL_TEST_CASE])
def test_create_zsl_interface(substrate, film, gap, vacuum, max_area, expected_interface):
    interface = create_zsl_interface(
        substrate_crystal=substrate.bulk_config,
        film_crystal=film.bulk_config,
        substrate_miller_indices=substrate.miller_indices,
        film_miller_indices=film.miller_indices,
        substrate_number_of_layers=substrate.number_of_layers,
        film_number_of_layers=film.number_of_layers,
        substrate_termination_formula=None,
        film_termination_formula=None,
        gap=gap,
        vacuum=vacuum,
        xy_shift=[0, 0],
        max_area=max_area,
        max_area_ratio_tol=MAX_AREA_RATIO_TOL,
        max_length_tol=MAX_LENGTH_TOL,
        max_angle_tol=MAX_ANGLE_TOL,
        use_conventional_cell=True,
        reduce_result_cell=False,
        reduce_result_cell_to_primitive=True,
    )
    interface.metadata.build = []
    expected_interface = get_platform_specific_value(expected_interface)
    assert_two_entities_deep_almost_equal(interface, expected_interface)


@pytest.mark.parametrize("substrate, film, gap, vacuum, max_area, expected_interface", [GRAPHENE_NICKEL_TEST_CASE])
def test_create_zsl_interface_between_slabs(substrate, film, gap, vacuum, max_area, expected_interface):
    substrate_slab_config = SlabConfiguration.from_parameters(
        material_or_dict=substrate.bulk_config,
        miller_indices=substrate.miller_indices,
        number_of_layers=substrate.number_of_layers,
        vacuum=0.0,
        termination_formula=None,
        use_conventional_cell=True,
    )
    film_slab_config = SlabConfiguration.from_parameters(
        material_or_dict=film.bulk_config,
        miller_indices=film.miller_indices,
        number_of_layers=film.number_of_layers,
        vacuum=0.0,
        termination_formula=None,
        use_conventional_cell=True,
    )

    substrate_slab = SlabBuilder().get_material(substrate_slab_config)
    film_slab = SlabBuilder().get_material(film_slab_config)

    interface = create_zsl_interface_between_slabs(
        substrate_slab=substrate_slab,
        film_slab=film_slab,
        gap=gap,
        vacuum=vacuum,
        xy_shift=[0, 0],
        max_area=max_area,
        max_area_ratio_tol=MAX_AREA_RATIO_TOL,
        max_length_tol=MAX_LENGTH_TOL,
        max_angle_tol=MAX_ANGLE_TOL,
        reduce_result_cell=False,
        reduce_result_cell_to_primitive=True,
    )
    expected_interface = get_platform_specific_value(expected_interface)
    assert_two_entities_deep_almost_equal(interface, expected_interface)
