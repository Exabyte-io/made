import functools
import types
import numpy as np
from typing import Union, List
from enum import Enum
from pymatgen.analysis.interfaces.coherent_interfaces import Interface
from pymatgen.core.structure import Structure


def patch_interface_with_mean_abs_strain(target: Interface, tolerance: float = 10e-6):
    def get_mean_abs_strain(target):
        return target.interface_properties[StrainModes.mean_abs_strain]

    target.get_mean_abs_strain = types.MethodType(get_mean_abs_strain, target)
    target.interface_properties[StrainModes.mean_abs_strain] = (
        round(np.mean(np.abs(target.interface_properties["strain"])) / tolerance) * tolerance
    )
    return target


class StrainModes(Enum):
    strain = "strain"
    von_mises_strain = "von_mises_strain"
    mean_abs_strain = "mean_abs_strain"


class InterfaceDataHolder(object):
    """
    A class to hold data for interfaces generated by pymatgen.
    Structures are stored in a dictionary with the termination as the key.
    Example data structure:
        {
            "('C_P6/mmm_2', 'Si_R-3m_1')": [
                { ...interface for ('C_P6/mmm_2', 'Si_R-3m_1') at index 0...},
                { ...interface for ('C_P6/mmm_2', 'Si_R-3m_1') at index 1...},
                ...
            ],
            "<termination at index 1>": [
                { ...interface for 'termination at index 1' at index 0...},
                { ...interface for 'termination at index 1' at index 1...},
                ...
            ]
        }
    """

    def __init__(self) -> None:
        self.data = {}
        self.terminations = []

    def add_termination(self, termination: str):
        if termination not in self.terminations:
            self.terminations.append(termination)
            self.set_interfaces_for_termination(termination, [])

    def add_interfaces_for_termination(self, termination: str, interfaces: Union[List[Interface], Interface]):
        self.add_termination(termination)
        self.set_interfaces_for_termination(termination, self.get_interfaces_for_termination(termination) + interfaces)

    def add_data_entries(
        self, entries=[], sort_interfaces_for_all_terminations_by_strain_and_size=True, remove_duplicates=True
    ):
        if isinstance(entries, Interface):
            entries = [entries]
        all_terminations = [e.interface_properties["termination"] for e in entries]
        unique_terminations = list(set(all_terminations))
        for termination in unique_terminations:
            entries_for_termination = [
                entry for entry in entries if entry.interface_properties["termination"] == termination
            ]
            self.add_interfaces_for_termination(termination, entries_for_termination)
        if sort_interfaces_for_all_terminations_by_strain_and_size:
            self.sort_interfaces_for_all_terminations_by_strain_and_size()
        if remove_duplicates:
            self.remove_duplicate_interfaces()

    def set_interfaces_for_termination(self, termination, interfaces):
        self.data[termination] = interfaces

    def get_interfaces_for_termination(self, termination):
        if isinstance(termination, int):
            termination = self.terminations[termination]
        return self.data.get(termination, [])

    def remove_duplicate_interfaces(self, strain_mode=StrainModes.mean_abs_strain):
        for termination in self.terminations:
            self.remove_duplicate_interfaces_for_termination(termination, strain_mode)

    def remove_duplicate_interfaces_for_termination(self, termination, strain_mode=StrainModes.mean_abs_strain):
        def are_interfaces_duplicate(interface1, interface2):
            return interface1.num_sites == interface2.num_sites and np.allclose(
                interface1.interface_properties[strain_mode], interface2.interface_properties[strain_mode]
            )

        sorted_interfaces = self.get_interfaces_for_termination_sorted_by_size(termination)
        filtered_interfaces = [sorted_interfaces[0]] if sorted_interfaces else []

        for interface in sorted_interfaces[1:]:
            if not any(
                are_interfaces_duplicate(interface, unique_interface) for unique_interface in filtered_interfaces
            ):
                filtered_interfaces.append(interface)

        self.set_interfaces_for_termination(termination, filtered_interfaces)

    def get_interfaces_for_termination_sorted_by_strain(self, termination, strain_mode=StrainModes.mean_abs_strain):
        return sorted(
            self.get_interfaces_for_termination(termination),
            key=lambda x: np.mean(np.abs(x.interface_properties[strain_mode])),
        )

    def get_interfaces_for_termination_sorted_by_size(self, termination):
        return sorted(
            self.get_interfaces_for_termination(termination),
            key=lambda x: x.num_sites,
        )

    def get_interfaces_for_termination_sorted_by_strain_and_size(
        self, termination, strain_mode=StrainModes.mean_abs_strain
    ):
        return sorted(
            self.get_interfaces_for_termination_sorted_by_strain(termination, strain_mode),
            key=lambda x: x.num_sites,
        )

    def sort_interfaces_for_all_terminations_by_strain_and_size(self):
        for termination in self.terminations:
            self.set_interfaces_for_termination(
                termination, self.get_interfaces_for_termination_sorted_by_strain_and_size(termination)
            )

    def get_all_interfaces(self):
        return functools.reduce(lambda a, b: a + b, self.data.values())

    def get_mean_abs_strain_for_interface(self, interface: Interface, tolerance: float = 10e-6) -> float:
        return round(np.mean(np.abs(interface.interface_properties["strain"])) / tolerance) * tolerance
