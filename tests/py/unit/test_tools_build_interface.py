import platform

from mat3ra.made.tools.build.interface import InterfaceConfiguration, ZSLInterfaceBuilder, ZSLStrainMatchingParameters

from .fixtures import INTERFACE_TERMINATION_PAIR, LAYER_CONFIGURATION, SUBSTRATE_CONFIGURATION

MAX_AREA = 100
# pymatgen `2023.6.23` supporting py3.8 returns 1 interface instead of 2
EXPECTED_NUMBER_OF_INTERFACES = 1 if platform.python_version().startswith("3.8") else 1
interface_configuration = InterfaceConfiguration(
    substrate_configuration=SUBSTRATE_CONFIGURATION,
    film_configuration=LAYER_CONFIGURATION,
    termination_pair=INTERFACE_TERMINATION_PAIR,
    distance_z=3.0,
    vacuum=5.0,
)

strain_matching_parameters = ZSLStrainMatchingParameters(max_area=MAX_AREA)


def test_create_interfaces():
    interfaces = ZSLInterfaceBuilder(strain_matching_parameters=strain_matching_parameters).get_materials(
        interface_configuration
    )
    assert len(interfaces) == EXPECTED_NUMBER_OF_INTERFACES
