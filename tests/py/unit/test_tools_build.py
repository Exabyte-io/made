import platform

from mat3ra.made.tools.build import create_interfaces
from mat3ra.made.tools.build.interface import InterfaceSettings

from .fixtures import LAYER_MATERIAL, SUBSTRATE_MATERIAL

MAX_AREA = 200
# pymatgen `2023.6.23` supporting py3.8 returns 1 interface instead of 2
EXPECTED_NUMBER_OF_INTERFACES = 1 if platform.python_version().startswith("3.8") else 2
settings = InterfaceSettings(
    USE_CONVENTIONAL_CELL=True,
    INTERFACE_PARAMETERS={"DISTANCE_Z": 3.0, "MAX_AREA": MAX_AREA},
    ZSL_PARAMETERS={"MAX_AREA": MAX_AREA, "MAX_AREA_TOL": 0.09, "MAX_LENGTH_TOL": 0.03, "MAX_ANGLE_TOL": 0.01},
    SUBSTRATE_PARAMETERS={"MILLER_INDICES": (1, 1, 1), "THICKNESS": 3},
    LAYER_PARAMETERS={"MILLER_INDICES": (0, 0, 1), "THICKNESS": 1},
)


def test_create_interfaces():
    interfaces = create_interfaces(substrate=SUBSTRATE_MATERIAL, layer=LAYER_MATERIAL, settings=settings)
    assert len(interfaces.get_interfaces_for_termination(0)) == EXPECTED_NUMBER_OF_INTERFACES
