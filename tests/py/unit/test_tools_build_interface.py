import platform

import pytest
from mat3ra.made.material import Material
from mat3ra.made.tools.build.interface import (
    ZSLStrainMatchingInterfaceBuilder,
    ZSLStrainMatchingInterfaceBuilderParameters,
    ZSLStrainMatchingParameters,
    create_interfaces,
)
from mat3ra.made.tools.build.interface.builders import (
    CommensurateLatticeTwistedInterfaceBuilder,
    CommensurateLatticeTwistedInterfaceBuilderParameters,
    NanoRibbonTwistedInterfaceBuilder,
    NanoRibbonTwistedInterfaceConfiguration,
    TwistedInterfaceConfiguration,
)
from mat3ra.utils import assertion as assertion_utils

from .fixtures.monolayer import GRAPHENE

MAX_AREA = 100
# pymatgen `2023.6.23` supporting py3.8 returns 1 interface instead of 2
EXPECTED_NUMBER_OF_INTERFACES = 1 if platform.python_version().startswith("3.8") else 1
interface_configuration = None

zsl_strain_matching_parameters = ZSLStrainMatchingParameters(max_area=MAX_AREA)
build_parameters = ZSLStrainMatchingInterfaceBuilderParameters(
    strain_matching_parameters=zsl_strain_matching_parameters
)
matched_interfaces_builder = ZSLStrainMatchingInterfaceBuilder(build_parameters=build_parameters)


@pytest.mark.skip(reason="Fixtures are commented out. To be fixed in epic-7623")
def test_create_interfaces():
    interfaces = create_interfaces(matched_interfaces_builder, interface_configuration)

    assert len(interfaces) == EXPECTED_NUMBER_OF_INTERFACES
    assert interfaces[0].name is not None


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
