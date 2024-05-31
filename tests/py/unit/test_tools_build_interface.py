import platform

from mat3ra.made.tools.build.interface import (
    InterfaceConfiguration,
    ZSLStrainMatchingInterfaceBuilder,
    ZSLStrainMatchingInterfaceBuilderParameters,
    ZSLStrainMatchingParameters,
)

from .fixtures import INTERFACE_TERMINATION_PAIR, LAYER_CONFIGURATION, SUBSTRATE_CONFIGURATION

MAX_AREA = 100
# pymatgen `2023.6.23` supporting py3.8 returns 1 interface instead of 2
EXPECTED_NUMBER_OF_INTERFACES = 1 if platform.python_version().startswith("3.8") else 1
interface_configuration = InterfaceConfiguration(
    film_configuration=LAYER_CONFIGURATION,
    substrate_configuration=SUBSTRATE_CONFIGURATION,
    film_termination=INTERFACE_TERMINATION_PAIR[0],
    substrate_termination=INTERFACE_TERMINATION_PAIR[1],
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
    interfaces = matched_interfaces_builder.get_materials(configuration=interface_configuration)

    assert len(interfaces) == EXPECTED_NUMBER_OF_INTERFACES
