from typing import List

import numpy as np

from mat3ra.made.material import Material
from .builders import InterfaceBuilder
from .configuration import InterfaceConfiguration
from .enums import StrainModes
from ..vacuum.configuration import VacuumConfiguration
from ...convert import PymatgenInterface, INTERFACE_LABELS_MAP
from ...modify import filter_by_label


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


def _build_interface_from_stack(
    selected_config,
    vacuum: float,
    gaps=None,
    xy_shift=None,
    direction=None,
) -> Material:
    vacuum_configuration = VacuumConfiguration(size=vacuum)
    stack_components = [
        selected_config.substrate_configuration,
        selected_config.film_configuration,
        vacuum_configuration,
    ]
    interface_config_kwargs = {
        "stack_components": stack_components,
    }
    if gaps is not None:
        interface_config_kwargs["gaps"] = gaps
    if xy_shift is not None:
        interface_config_kwargs["xy_shift"] = xy_shift
    if direction is not None:
        interface_config_kwargs["direction"] = direction
    interface_config = InterfaceConfiguration(**interface_config_kwargs)
    builder = InterfaceBuilder()
    return builder.get_material(interface_config)
