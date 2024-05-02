from mat3ra.made.tools.build import create_interfaces
from mat3ra.made.tools.build.interface import InterfaceSettings
import os
import json

dir_path = os.path.dirname(os.path.realpath(__file__))

substrate_path = os.path.join(dir_path, "../../fixtures/Ni-hex.json")
layer_path = os.path.join(dir_path, "../../fixtures/Graphene.json")

with open(substrate_path) as file:
    substrate_material = json.load(file)

with open(layer_path) as file:
    layer_material = json.load(file)

settings = InterfaceSettings(
    USE_CONVENTIONAL_CELL=True,
    INTERFACE_PARAMETERS={"MAX_AREA": 50, "DISTANCE_Z": 1.0},
    LAYER_PARAMETERS={"THICKNESS": 1},
    SUBSTRATE_PARAMETERS={"THICKNESS": 1},
)


def test_create_interfaces():
    interfaces = create_interfaces(substrate_material, layer_material, settings)
    assert len(interfaces.get_interfaces_for_termination(0)) == 1
