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

from . import BaseBuilder
from ..convert import to_ase, from_ase, to_pymatgen, from_pymatgen
from .slab import BaseSlabConfiguration, SlabConfiguration
from ...material import Material

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
        # TODO: implemetn utils to remove vacuum
        return self.get_interface()

    def get_interface(self) -> Material:
        return self.__builder.get_material(self)


class InterfaceBuilder(BaseBuilder):
    _ConfigurationType = InterfaceConfiguration


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


class ZSLStrainMatchingInterfaceBuilder(StrainMatchingInterfaceBuilder):
    """
    Creates matching interface between substrate and film using the ZSL algorithm.
    """

    _BuildParametersType = ZSLStrainMatchingInterfaceBuilderParameters
    _GeneratedItemType = PymatgenInterface

    def _generate(self, configuration: InterfaceConfiguration) -> List[PymatgenInterface]:
        generator = ZSLGenerator(**self.build_parameters.dict())
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
        # TODO: sort by number of atoms
        return sorted(items, key=lambda x: np.mean(np.abs(x.interface_properties[StrainModes.mean_abs_strain])))

    def _post_process(self, items: List[_GeneratedItemType], post_process_parameters=None) -> List[Material]:
        # TODO: add metadata, and change name
        return [Material(from_pymatgen(interface)) for interface in items]


def interface_patch_with_mean_abs_strain(target: PymatgenInterface, tolerance: float = 10e-6):
    def get_mean_abs_strain(target):
        return target.interface_properties[StrainModes.mean_abs_strain]

    target.get_mean_abs_strain = types.MethodType(get_mean_abs_strain, target)
    target.interface_properties[StrainModes.mean_abs_strain] = (
        round(np.mean(np.abs(target.interface_properties["strain"])) / tolerance) * tolerance
    )
    return target
