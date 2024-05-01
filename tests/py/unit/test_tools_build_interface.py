from mat3ra.made.tools.build.interface import InterfaceDataHolder
from mat3ra.made.tools.build.interface import patch_interface_with_mean_abs_strain, StrainModes

from .utils import INTERFACE_STRUCTURE

INTERFACE_STRUCTURE.interface_properties["termination"] = ("Si_termination", "Cu_termination")
INTERFACE_STRUCTURE.interface_properties[StrainModes.strain] = 0.1
INTERFACE_STRUCTURE = patch_interface_with_mean_abs_strain(INTERFACE_STRUCTURE)


def test_add_data_entries():
    interfaces_data = InterfaceDataHolder()
    interfaces_data.add_data_entries(INTERFACE_STRUCTURE)
    assert len(interfaces_data.get_interfaces_for_termination(0)) == 1


def test_get_interfaces_for_termination():
    interfaces_data = InterfaceDataHolder()
    interfaces_data.add_data_entries([INTERFACE_STRUCTURE])
    assert interfaces_data.get_interfaces_for_termination(0)[0] == INTERFACE_STRUCTURE
