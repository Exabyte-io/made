import functools
import types
import numpy as np
from typing import Union, List, Tuple, Dict, Optional
from enum import Enum
from pydantic import BaseModel
from mat3ra.utils import array as array_utils
from pymatgen.analysis.interfaces.coherent_interfaces import CoherentInterfaceBuilder, ZSLGenerator, Interface

from ..modify import wrap_to_unit_cell
from ..convert import convert_atoms_or_structure_to_material, to_ase, from_ase
from .slab import BaseSlabConfiguration, SlabConfiguration

TerminationPair = Tuple[str, str]
InterfacesType = List[Interface]
InterfacesDataType = Dict[Tuple, List[Interface]]


class StrainModes(str, Enum):
    strain = "strain"
    von_mises_strain = "von_mises_strain"
    mean_abs_strain = "mean_abs_strain"


class ZSLStrainMatchingParameters(BaseModel):
    max_area: float = 50.0
    max_area_tol: float = 0.09
    max_length_tol: float = 0.03
    max_angle_tol: float = 0.01


class InterfaceConfiguration(BaseSlabConfiguration):
    def __init__(
        self,
        substrate_configuration: SlabConfiguration,
        film_configuration: SlabConfiguration,
        termination_pair: TerminationPair,
        distance_z: float = 3.0,
        vacuum: float = 5.0,
        shift_x: float = 0.0,
        shift_y: float = 0.0,
    ):
        super().__init__()
        self.substrate_configuration = substrate_configuration
        self.film_configuration = film_configuration
        self.termination_pair = termination_pair
        self.vacuum: float = vacuum
        self.distance_z: float = distance_z
        self.shift_x: float = shift_x
        self.shift_y: float = shift_y

    @property
    def bulk(self):
        # TODO: Return the bulk structure of the interface configuration, with no vacuum
        # Slab with vacuum=0 and thickness=1 ?
        return self.get_material()

    @property
    def miller_indices(self):
        # Return the bulk structure of the interface configuration, with no vacuum
        # Thickness should be one
        return self.__miller_indices

    def get_material(self):
        substrate_slab = self.substrate_configuration.get_material(termination=self.termination_pair[0])
        film_slab = self.film_configuration.get_material(termination=self.termination_pair[1])

        substrate_slab_ase = to_ase(substrate_slab)
        film_slab_ase = to_ase(film_slab)

        # Calculate z-shift
        max_z_substrate = max(substrate_slab_ase.positions[:, 2])
        min_z_film = min(film_slab_ase.positions[:, 2])
        z_shift = max_z_substrate - min_z_film + self.distance_z

        # Apply z-shift to the film
        film_slab_ase.translate((0, 0, z_shift))

        # TODO: add x,y shift: in cartesian and crystal coordinates

        # Combine the substrate and film into one Atoms object
        interface_ase = substrate_slab_ase + film_slab_ase

        new_cell_height = max(interface_ase.positions[:, 2]) + self.vacuum
        interface_ase.cell[2, 2] = new_cell_height

        return from_ase(interface_ase)

    @property
    def interface_properties(self):
        # Strain
        pass

    def termination_pairs(self):
        # TODO
        pass


class InterfaceConfigurationStrainMatcher:
    def __init__(
        self,
        interface_configuration: InterfaceConfiguration,
        strain_matching_parameters: Optional[ZSLStrainMatchingParameters],
        termination_pair: TerminationPair,
    ):
        self.termination_pair = termination_pair
        self.interface_configuration = interface_configuration
        self.strain_matching_parameters = strain_matching_parameters

    @property
    def generator(self):
        if isinstance(self.strain_matching_parameters, ZSLStrainMatchingParameters):
            zsl_parameters = ZSLStrainMatchingParameters(**self.strain_matching_parameters.__dict__)
            generator: ZSLGenerator = ZSLGenerator(**zsl_parameters)
            return CoherentInterfaceBuilder(
                substrate_structure=self.interface_configuration.substrate_configuration.bulk,
                film_structure=self.interface_configuration.film_configuration.bulk,
                substrate_miller=self.interface_configuration.substrate_configuration.miller_indices,
                film_miller=self.interface_configuration.film_configuration.miller_indices,
                zslgen=generator,
            )
        elif self.strain_matching_parameters is None:

            class ScaleInterfaceBuilder:
                def __init__(self, interface_configuration: InterfaceConfiguration):
                    super().__init__(
                        interface_configuration=interface_configuration,
                        film_structure=interface_configuration.film_configuration.bulk,
                        substrate_miller=interface_configuration.substrate_configuration.miller_indices,
                        film_miller=interface_configuration.film_configuration.miller_indices,
                    )

                def get_interfaces(self):
                    interface = ...
                    return [interface]

            return ScaleInterfaceBuilder(self.interface_configuration)
        else:
            raise ValueError(f"Unsupported strain matching algorithm: {self.strain_matching_algorithm}")

    # def builder(self):
    #     return CoherentInterfaceBuilder(
    #         substrate_structure=self.interface_configuration.substrate_configuration.bulk,
    #         film_structure=self.interface_configuration.film_configuration.bulk,
    #         substrate_miller=self.interface_configuration.substrate_configuration.miller_indices,
    #         film_miller=self.interface_configuration.film_configuration.miller_indices,
    #         zslgen=self.generator,
    #     )

    def get_interface_data_holder(self):
        builder = self.generator
        interfaces_data = InterfaceDataHolder()

        all_interfaces_for_termination = builder.get_interfaces(
            termination=self.termination_pair,
            gap=self.interface_configuration.distance_z,
            film_thickness=self.interface_configuration.film_configuration.thickness,
            substrate_thickness=self.interface_configuration.substrate_configuration.thickness,
            in_layers=True,
        )

        all_interfaces_for_termination_patched_wrapped = list(
            map(
                lambda i: wrap_to_unit_cell(interface_patch_with_mean_abs_strain(i)),
                all_interfaces_for_termination,
            )
        )

        interfaces_data.add_data_entries(
            all_interfaces_for_termination_patched_wrapped,
            sort_interfaces_by_strain_and_size=True,
            remove_duplicates=True,
        )


def interface_patch_with_mean_abs_strain(target: Interface, tolerance: float = 10e-6):
    def get_mean_abs_strain(target):
        return target.interface_properties[StrainModes.mean_abs_strain]

    target.get_mean_abs_strain = types.MethodType(get_mean_abs_strain, target)
    target.interface_properties[StrainModes.mean_abs_strain] = (
        round(np.mean(np.abs(target.interface_properties["strain"])) / tolerance) * tolerance
    )
    return target


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

    def __init__(self, entries: Union[InterfacesType, None] = None) -> None:
        if entries is None:
            entries = []
        self.data: InterfacesDataType = {}
        self.terminations: List[TerminationPair] = []
        self.add_data_entries(entries)

    def __str__(self):
        terminations_list = f"There are {len(self.terminations)} terminations:" + ", ".join(
            f"\n{idx}: ({a}, {b})" for idx, (a, b) in enumerate(self.terminations)
        )
        interfaces_list = "\n".join(
            [
                f"There are {len(self.data[termination])} interfaces for termination {termination}:\n{idx}: "
                + f"{self.data[termination]}"
                for idx, termination in enumerate(self.terminations)
            ]
        )
        return f"{terminations_list}\n{interfaces_list}"

    def add_termination(self, termination: Tuple[str, str]):
        if termination not in self.terminations:
            self.terminations.append(termination)
            self.set_interfaces_for_termination(termination, [])

    def add_interfaces_for_termination(
        self, termination: TerminationPair, interfaces: Union[InterfacesType, Interface]
    ):
        self.add_termination(termination)
        self.set_interfaces_for_termination(termination, self.get_interfaces_for_termination(termination) + interfaces)

    def add_data_entries(
        self,
        entries: List[Interface] = [],
        sort_interfaces_by_strain_and_size: bool = True,
        remove_duplicates: bool = True,
    ):
        entries = array_utils.convert_to_array_if_not(entries)
        all_terminations = [e.interface_properties["termination"] for e in entries]
        unique_terminations = list(set(all_terminations))
        for termination in unique_terminations:
            entries_for_termination = [
                entry for entry in entries if entry.interface_properties["termination"] == termination
            ]
            self.add_interfaces_for_termination(termination, entries_for_termination)
        if sort_interfaces_by_strain_and_size:
            self.sort_interfaces_for_all_terminations_by_strain_and_size()
        if remove_duplicates:
            self.remove_duplicate_interfaces()

    def set_interfaces_for_termination(self, termination: TerminationPair, interfaces: List[Interface]):
        self.data[termination] = interfaces

    def get_termination(self, termination: Union[int, TerminationPair]) -> TerminationPair:
        if isinstance(termination, int):
            termination = self.terminations[termination]
        return termination

    def get_interfaces_for_termination_or_its_index(
        self, termination_or_its_index: Union[int, TerminationPair]
    ) -> List[Interface]:
        termination = self.get_termination(termination_or_its_index)
        return self.data[termination]

    def get_interfaces_for_termination(
        self,
        termination_or_its_index: Union[int, TerminationPair],
        slice_or_index_or_indices: Union[int, slice, List[int], None] = None,
    ) -> List[Interface]:
        interfaces = self.get_interfaces_for_termination_or_its_index(termination_or_its_index)
        return array_utils.filter_by_slice_or_index_or_indices(interfaces, slice_or_index_or_indices)

    def remove_duplicate_interfaces(self, strain_mode: StrainModes = StrainModes.mean_abs_strain):
        for termination in self.terminations:
            self.remove_duplicate_interfaces_for_termination(termination, strain_mode)

    def remove_duplicate_interfaces_for_termination(
        self, termination, strain_mode: StrainModes = StrainModes.mean_abs_strain
    ):
        def are_interfaces_duplicate(interface1: Interface, interface2: Interface):
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

    def get_interfaces_for_termination_sorted_by_strain(
        self, termination: Union[int, TerminationPair], strain_mode: StrainModes = StrainModes.mean_abs_strain
    ) -> List[Interface]:
        return sorted(
            self.get_interfaces_for_termination(termination),
            key=lambda x: np.mean(np.abs(x.interface_properties[strain_mode])),
        )

    def get_interfaces_for_termination_sorted_by_size(
        self, termination: Union[int, TerminationPair]
    ) -> List[Interface]:
        return sorted(
            self.get_interfaces_for_termination(termination),
            key=lambda x: x.num_sites,
        )

    def get_interfaces_for_termination_sorted_by_strain_and_size(
        self, termination: Union[int, TerminationPair], strain_mode: StrainModes = StrainModes.mean_abs_strain
    ) -> List[Interface]:
        return sorted(
            self.get_interfaces_for_termination_sorted_by_strain(termination, strain_mode),
            key=lambda x: x.num_sites,
        )

    def sort_interfaces_for_all_terminations_by_strain_and_size(self):
        for termination in self.terminations:
            self.set_interfaces_for_termination(
                termination, self.get_interfaces_for_termination_sorted_by_strain_and_size(termination)
            )

    def get_all_interfaces(self) -> List[Interface]:
        return functools.reduce(lambda a, b: a + b, self.data.values())

    def get_interfaces_as_materials(
        self, termination: Union[int, TerminationPair], slice_range_or_index: Union[int, slice]
    ) -> List[Interface]:
        return list(
            map(
                convert_atoms_or_structure_to_material,
                self.get_interfaces_for_termination(termination, slice_range_or_index),
            )
        )
