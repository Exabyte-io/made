import types
from typing import List
import numpy as np
from mat3ra.made.material import Material

from ...modify import filter_by_label
from ...convert import PymatgenInterface, INTERFACE_LABELS_MAP
from .enums import StrainModes


def interface_patch_with_mean_abs_strain(target: PymatgenInterface, tolerance: float = 10e-6):
    def get_mean_abs_strain(target):
        return target.interface_properties[StrainModes.mean_abs_strain]

    target.get_mean_abs_strain = types.MethodType(get_mean_abs_strain, target)
    target.interface_properties[StrainModes.mean_abs_strain] = (
        round(np.mean(np.abs(target.interface_properties["strain"])) / tolerance) * tolerance
    )
    return target


def remove_duplicate_interfaces(
    interfaces: List[PymatgenInterface], strain_mode: StrainModes = StrainModes.mean_abs_strain
):
    def are_interfaces_duplicate(interface1: PymatgenInterface, interface2: PymatgenInterface):
        are_sizes_equivalent = interface1.num_sites == interface2.num_sites and np.allclose(
            interface1.interface_properties[strain_mode], interface2.interface_properties[strain_mode]
        )
        are_strains_equivalent = np.allclose(
            interface1.interface_properties[strain_mode], interface2.interface_properties[strain_mode]
        )
        return are_sizes_equivalent and are_strains_equivalent

    filtered_interfaces = [interfaces[0]] if interfaces else []

    for interface in interfaces[1:]:
        if not any(are_interfaces_duplicate(interface, unique_interface) for unique_interface in filtered_interfaces):
            filtered_interfaces.append(interface)
    return filtered_interfaces


def get_slab(interface: Material, part: str = "film"):
    try:
        return filter_by_label(interface, INTERFACE_LABELS_MAP[part])
    except ValueError:
        raise ValueError(f"Material does not contain label for {part}.")
