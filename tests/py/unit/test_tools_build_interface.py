import platform

from mat3ra.made.material import Material
from mat3ra.made.tools.build.interface import (
    InterfaceConfiguration,
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
from unit.fixtures.generated.fixtures import (
    FILM_CONFIGURATION,
    INTERFACE_NAME,
    INTERFACE_TERMINATION_PAIR,
    SUBSTRATE_CONFIGURATION,
)

from .fixtures.monolayer import GRAPHENE

MAX_AREA = 100
# pymatgen `2023.6.23` supporting py3.8 returns 1 interface instead of 2
EXPECTED_NUMBER_OF_INTERFACES = 1 if platform.python_version().startswith("3.8") else 1
interface_configuration = InterfaceConfiguration(
    film_configuration=FILM_CONFIGURATION,
    substrate_configuration=SUBSTRATE_CONFIGURATION,
    film_termination=INTERFACE_TERMINATION_PAIR.film_termination,
    substrate_termination=INTERFACE_TERMINATION_PAIR.substrate_termination,
    distance_z=3.0,
    vacuum=5.0,
)

zsl_strain_matching_parameters = ZSLStrainMatchingParameters(max_area=MAX_AREA)
build_parameters = ZSLStrainMatchingInterfaceBuilderParameters(
    strain_matching_parameters=zsl_strain_matching_parameters
)
matched_interfaces_builder = ZSLStrainMatchingInterfaceBuilder(build_parameters=build_parameters)


def test_create_interfaces():
    interfaces = create_interfaces(matched_interfaces_builder, interface_configuration)

    assert len(interfaces) == EXPECTED_NUMBER_OF_INTERFACES
    assert interfaces[0].name == INTERFACE_NAME


def test_create_twisted_nanoribbon_interface():
    film = Material.create(GRAPHENE)
    configuration = NanoRibbonTwistedInterfaceConfiguration(
        film=film,
        substrate=film,
        twist_angle=60,
        distance_z=3.0,
        ribbon_width=3,
        ribbon_length=5,
        vacuum_x=2.0,
        vacuum_y=2.0,
    )

    builder = NanoRibbonTwistedInterfaceBuilder()
    interface = builder.get_material(configuration)

    expected_cell_vectors = [[15.102811, 0.0, 0.0], [0.0, 16.108175208, 0.0], [0.0, 0.0, 20.0]]
    expected_coordinate = [0.704207885, 0.522108183, 0.65]
    assertion_utils.assert_deep_almost_equal(expected_cell_vectors, interface.lattice.vector_arrays)
    assertion_utils.assert_deep_almost_equal(expected_coordinate, interface.basis.coordinates.values[42])


def test_create_commensurate_supercell_twisted_interface():
    film = Material.create(GRAPHENE)
    substrate = Material.create(GRAPHENE)
    config = TwistedInterfaceConfiguration(film=film, substrate=substrate, twist_angle=13, distance_z=3.0)
    params = CommensurateLatticeTwistedInterfaceBuilderParameters(
        max_supercell_matrix_int=5, angle_tolerance=0.5, return_first_match=True
    )
    builder = CommensurateLatticeTwistedInterfaceBuilder(build_parameters=params)
    interfaces = builder.get_materials(config, post_process_parameters=config)
    assert len(interfaces) == 1
    interface = interfaces[0]
    expected_cell_vectors = [[10.754672133, 0.0, 0.0], [5.377336066500001, 9.313819276550575, 0.0], [0.0, 0.0, 20.0]]
    assertion_utils.assert_deep_almost_equal(expected_cell_vectors, interface.basis.cell.vector_arrays)
    expected_angle = 13.174
    assert interface.metadata["build"]["configuration"]["actual_twist_angle"] == expected_angle
