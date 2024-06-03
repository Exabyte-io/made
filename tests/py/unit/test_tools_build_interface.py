import platform

from mat3ra.made.tools.build.interface import (
    InterfaceConfiguration,
    ZSLStrainMatchingInterfaceBuilder,
    ZSLStrainMatchingInterfaceBuilderParameters,
    ZSLStrainMatchingParameters,
    create_interfaces,
)

from .fixtures import FILM_CONFIGURATION, INTERFACE_NAME, INTERFACE_TERMINATION_PAIR, SUBSTRATE_CONFIGURATION

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
matched_interfaces_builder = ZSLStrainMatchingInterfaceBuilder(
    build_parameters=ZSLStrainMatchingInterfaceBuilderParameters(
        strain_matching_parameters=zsl_strain_matching_parameters
    )
)


def test_create_interfaces():
    interfaces = create_interfaces(matched_interfaces_builder, interface_configuration)

    assert len(interfaces) == EXPECTED_NUMBER_OF_INTERFACES
    assert interfaces[0].name == INTERFACE_NAME
