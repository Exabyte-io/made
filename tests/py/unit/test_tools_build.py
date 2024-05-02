import json
import os

from mat3ra.made.material import Material
from mat3ra.made.tools.build import create_interfaces
from mat3ra.made.tools.build.interface import InterfaceSettings

dir_path = os.path.dirname(os.path.realpath(__file__))

substrate_path = os.path.join(dir_path, "../../fixtures/Ni-hex.json")
layer_path = os.path.join(dir_path, "../../fixtures/Graphene.json")

with open(substrate_path) as file:
    substrate_material = Material(json.load(file))

with open(layer_path) as file:
    layer_material = Material(json.load(file))


MAX_AREA = 200
EXPECTED_NUMBER_OF_INTERFACES = 8
settings = InterfaceSettings(
    USE_CONVENTIONAL_CELL=True,
    INTERFACE_PARAMETERS={"DISTANCE_Z": 3.0, "MAX_AREA": MAX_AREA},
    ZSL_PARAMETERS={"MAX_AREA": MAX_AREA, "MAX_AREA_TOL": 0.09, "MAX_LENGTH_TOL": 0.03, "MAX_ANGLE_TOL": 0.01},
    SUBSTRATE_PARAMETERS={"MILLER_INDICES": (1, 1, 1), "THICKNESS": 3},
    LAYER_PARAMETERS={"MILLER_INDICES": (0, 0, 1), "THICKNESS": 1},
)


def test_create_interfaces():
    interfaces = create_interfaces(substrate_material, layer_material, settings)
    assert len(interfaces.get_interfaces_for_termination(0)) == EXPECTED_NUMBER_OF_INTERFACES
