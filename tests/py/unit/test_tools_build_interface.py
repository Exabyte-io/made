from mat3ra.made.tools.build.interface import InterfaceDataHolder

from .fixtures import INTERFACE_STRUCTURE, INTERFACE_TERMINATION


def test_add_data_entries():
    interfaces_data = InterfaceDataHolder()
    interfaces_data.add_data_entries(INTERFACE_STRUCTURE)
    assert len(interfaces_data.get_interfaces_for_termination(0)) == 1
    assert len(interfaces_data.get_interfaces_for_termination(INTERFACE_TERMINATION)) == 1


def test_get_interfaces_for_termination():
    interfaces_data = InterfaceDataHolder()
    interfaces_data.add_data_entries([INTERFACE_STRUCTURE])
    assert interfaces_data.get_interfaces_for_termination(0)[0] == INTERFACE_STRUCTURE


def test_remove_duplicate_interfaces():
    interfaces_data = InterfaceDataHolder()
    interfaces_data.add_data_entries([INTERFACE_STRUCTURE, INTERFACE_STRUCTURE], remove_duplicates=False)
    assert len(interfaces_data.get_interfaces_for_termination(INTERFACE_TERMINATION)) == 2
    interfaces_data.remove_duplicate_interfaces()
    assert len(interfaces_data.get_interfaces_for_termination(INTERFACE_TERMINATION)) == 1
