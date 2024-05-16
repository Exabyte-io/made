import platform

from mat3ra.made.tools.build import create_interfaces
from mat3ra.made.tools.build.interface import InterfaceBuilderSettings, SlabParameters, StrainMatchingParameters

from .fixtures import LAYER_MATERIAL, SUBSTRATE_MATERIAL

MAX_AREA = 200
# pymatgen `2023.6.23` supporting py3.8 returns 1 interface instead of 2
EXPECTED_NUMBER_OF_INTERFACES = 1 if platform.python_version().startswith("3.8") else 2
settings = InterfaceBuilderSettings(
    use_conventional_cell=True,
    substrate_parameters=SlabParameters(miller_indices=(1, 1, 1)),
    strain_matching_parameters=StrainMatchingParameters(max_area=MAX_AREA),
)


def test_create_interfaces():
    interfaces = create_interfaces(substrate=SUBSTRATE_MATERIAL, layer=LAYER_MATERIAL, settings=settings)
    assert len(interfaces.get_interfaces_for_termination(0)) == EXPECTED_NUMBER_OF_INTERFACES
