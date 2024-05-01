from mat3ra.made.tools.build.interface import InterfaceDataHolder
from .fixtures import INTERFACE_STRUCTURE


def test_add_data_entries():
    interfaces_data = InterfaceDataHolder()
    interfaces_data.add_data_entries(INTERFACE_STRUCTURE)
    assert len(interfaces_data.get_interfaces_for_termination(0)) == 1


def test_get_interfaces_for_termination():
    interfaces_data = InterfaceDataHolder()
    interfaces_data.add_data_entries([INTERFACE_STRUCTURE])
    assert interfaces_data.get_interfaces_for_termination(0)[0] == INTERFACE_STRUCTURE
