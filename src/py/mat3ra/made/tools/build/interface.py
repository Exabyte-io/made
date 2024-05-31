import types

import numpy as np
from typing import List, Tuple, Dict, Optional, Any
from enum import Enum
from pydantic import BaseModel
from pymatgen.analysis.interfaces.coherent_interfaces import (
    CoherentInterfaceBuilder,
    ZSLGenerator,
    Interface as PymatgenInterface,
)
from ase.build.tools import niggli_reduce
from ase import Atoms as ASEAtoms

from ...material import Material
from ..convert import to_ase, from_ase, to_pymatgen, from_pymatgen
from ..build import BaseBuilder
from .slab import BaseSlabConfiguration, SlabConfiguration

TerminationPair = Tuple[str, str]
InterfacesType = List[PymatgenInterface]
InterfacesDataType = Dict[Tuple, List[PymatgenInterface]]


class StrainModes(str, Enum):
    strain = "strain"
    von_mises_strain = "von_mises_strain"
    mean_abs_strain = "mean_abs_strain"


class SimpleInterfaceBuilderParameters(BaseModel):
    scale_film: bool = True


class StrainMatchingInterfaceBuilderParameters(BaseModel):
    strain_matching_parameters: Optional[Any] = None


class ZSLStrainMatchingParameters(BaseModel):
    max_area: float = 50.0
    max_area_ratio_tol: float = 0.09
    max_length_tol: float = 0.03
    max_angle_tol: float = 0.01


class ZSLStrainMatchingInterfaceBuilderParameters(StrainMatchingInterfaceBuilderParameters):
    strain_matching_parameters: ZSLStrainMatchingParameters


class CSLStrainMatchingParameters(BaseModel):
    pass


class CSLStrainMatchingInterfaceBuilderParameters(StrainMatchingInterfaceBuilderParameters):
    strain_matching_parameters: CSLStrainMatchingParameters


class InterfaceConfiguration(BaseSlabConfiguration):
    film_configuration: SlabConfiguration
    substrate_configuration: SlabConfiguration
    film_termination: str
    substrate_termination: str
    termination_pair: TerminationPair
    distance_z: float = 3.0
    vacuum: float = 5.0

    def __init__(
        self,
        film_configuration: SlabConfiguration,
        substrate_configuration: SlabConfiguration,
        film_termination: str,
        substrate_termination: str,
        distance_z: float = 3.0,
        vacuum: float = 5.0,
    ):
        super().__init__()
        self.film_configuration = film_configuration
        self.substrate_configuration = substrate_configuration
        self.film_termination = film_termination
        self.substrate_termination = substrate_termination
        self.termination_pair = (film_termination, substrate_termination)
        self.distance_z: float = distance_z
        self.vacuum: float = vacuum
        self.__builder = SimpleInterfaceBuilder(build_parameters=SimpleInterfaceBuilderParameters(scale_film=False))

    @property
    def bulk(self):
        # TODO: implement utils to remove vacuum
        return self.get_interface()

    def get_interface(self) -> Material:
        return self.__builder.get_material(self)


class InterfaceBuilder(BaseBuilder):
    _ConfigurationType: type(InterfaceConfiguration) = InterfaceConfiguration  # type: ignore

    def _finalize(self, materials: List[Material], configuration: _ConfigurationType) -> List[Material]:
        return [self._update_material_name(material, configuration) for material in materials]

    def _update_material_name(self, material: Material, configuration: InterfaceConfiguration) -> Material:
        film_formula = configuration.film_configuration.bulk["name"]
        substrate_formula = configuration.substrate_configuration.bulk["name"]
        film_miller_indices = "".join([str(i) for i in configuration.film_configuration.miller_indices])
        substrate_miller_indices = "".join([str(i) for i in configuration.substrate_configuration.miller_indices])
        new_name = f"{film_formula}({film_miller_indices})-{substrate_formula}({substrate_miller_indices}) Interface"
        material.name = new_name
        return material


class SimpleInterfaceBuilder(InterfaceBuilder):
    """
    Creates matching interface between substrate and film by straining the film to match the substrate.
    """

    _BuildParametersType = Optional[SimpleInterfaceBuilderParameters]
    _GeneratedItemType: ASEAtoms = ASEAtoms

    @staticmethod
    def __preprocess_slab_configuration(configuration: SlabConfiguration, termination: str):
        slab = configuration.get_slab(termination=termination)
        ase_slab = to_ase(slab)
        niggli_reduce(ase_slab)
        return ase_slab

    @staticmethod
    def __combine_two_slabs_ase(substrate_slab_ase: ASEAtoms, film_slab_ase: ASEAtoms, distance_z: float) -> ASEAtoms:
        max_z_substrate = max(substrate_slab_ase.positions[:, 2])
        min_z_film = min(film_slab_ase.positions[:, 2])
        shift_z = max_z_substrate - min_z_film + distance_z

        film_slab_ase.translate([0, 0, shift_z])

        return substrate_slab_ase + film_slab_ase

    @staticmethod
    def __add_vacuum_along_c_ase(interface_ase: ASEAtoms, vacuum: float) -> ASEAtoms:
        cell_c_with_vacuum = max(interface_ase.positions[:, 2]) + vacuum
        interface_ase.cell[2, 2] = cell_c_with_vacuum
        return interface_ase

    def _generate(self, configuration: InterfaceBuilder._ConfigurationType) -> List[_GeneratedItemType]:  # type: ignore
        film_slab_ase = self.__preprocess_slab_configuration(
            configuration.film_configuration, configuration.film_termination
        )
        substrate_slab_ase = self.__preprocess_slab_configuration(
            configuration.substrate_configuration, configuration.substrate_termination
        )

        if self.build_parameters.scale_film:
            film_slab_ase.set_cell(substrate_slab_ase.cell, scale_atoms=True)
            film_slab_ase.wrap()

        interface_ase = self.__combine_two_slabs_ase(substrate_slab_ase, film_slab_ase, configuration.distance_z)
        interface_ase_with_vacuum = self.__add_vacuum_along_c_ase(interface_ase, configuration.vacuum)

        return [interface_ase_with_vacuum]

    def _post_process(self, items: List[ASEAtoms], post_process_parameters=None) -> List[Material]:
        return [Material(from_ase(slab)) for slab in items]


class StrainMatchingInterfaceBuilder(InterfaceBuilder):
    _BuildParametersType = StrainMatchingInterfaceBuilderParameters

    def _update_material_name(self, material: Material, configuration: InterfaceConfiguration) -> Material:
        updated_material = super()._update_material_name(material, configuration)
        if StrainModes.mean_abs_strain in material.metadata:
            strain = material.metadata[StrainModes.mean_abs_strain]
            new_name = f"{updated_material.name}, Strain (mean_abs_strain={strain:.3f}%)"
            updated_material.name = new_name
        return material

    def _finalize(
        self, materials: List[Material], configuration: InterfaceBuilder._ConfigurationType
    ) -> List[Material]:
        return [self._update_material_name(material, configuration) for material in materials]


class ZSLStrainMatchingInterfaceBuilder(StrainMatchingInterfaceBuilder):
    """
    Creates matching interface between substrate and film using the ZSL algorithm.
    """

    _BuildParametersType: type(  # type: ignore
        ZSLStrainMatchingInterfaceBuilderParameters
    ) = ZSLStrainMatchingInterfaceBuilderParameters  # type: ignore
    _GeneratedItemType: PymatgenInterface = PymatgenInterface  # type: ignore

    def _generate(self, configuration: InterfaceConfiguration) -> List[PymatgenInterface]:
        generator = ZSLGenerator(**self.build_parameters.strain_matching_parameters.dict())
        builder = CoherentInterfaceBuilder(
            substrate_structure=to_pymatgen(configuration.substrate_configuration.bulk),
            film_structure=to_pymatgen(configuration.film_configuration.bulk),
            substrate_miller=configuration.substrate_configuration.miller_indices,
            film_miller=configuration.film_configuration.miller_indices,
            zslgen=generator,
        )

        interfaces = builder.get_interfaces(
            termination=configuration.termination_pair,
            gap=configuration.distance_z,
            film_thickness=configuration.film_configuration.thickness,
            substrate_thickness=configuration.substrate_configuration.thickness,
            in_layers=True,
        )

        return list([interface_patch_with_mean_abs_strain(interface) for interface in interfaces])

    def _sort(self, items: List[_GeneratedItemType]):
        sorted_by_num_sites = sorted(items, key=lambda x: x.num_sites)
        sorted_by_num_sites_and_strain = sorted(
            sorted_by_num_sites, key=lambda x: np.mean(x.interface_properties[StrainModes.mean_abs_strain])
        )
        unique_sorted_interfaces = remove_duplicate_interfaces(
            sorted_by_num_sites_and_strain, strain_mode=StrainModes.mean_abs_strain
        )
        return unique_sorted_interfaces

    def _post_process(self, items: List[_GeneratedItemType], post_process_parameters=None) -> List[Material]:
        materials = [Material(from_pymatgen(interface)) for interface in items]
        strains = [interface.interface_properties[StrainModes.mean_abs_strain] for interface in items]

        for material, strain in zip(materials, strains):
            material.metadata["mean_abs_strain"] = strain
        return materials


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
        size_the_same = interface1.num_sites == interface2.num_sites and np.allclose(
            interface1.interface_properties[strain_mode], interface2.interface_properties[strain_mode]
        )
        strain_the_same = np.allclose(
            interface1.interface_properties[strain_mode], interface2.interface_properties[strain_mode]
        )
        return size_the_same and strain_the_same

    filtered_interfaces = [interfaces[0]] if interfaces else []

    for interface in interfaces[1:]:
        if not any(are_interfaces_duplicate(interface, unique_interface) for unique_interface in filtered_interfaces):
            filtered_interfaces.append(interface)
    return filtered_interfaces
