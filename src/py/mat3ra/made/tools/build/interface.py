import numpy as np
from typing import List, Tuple, Dict, Optional, Union
from enum import Enum
from pydantic import BaseModel
from pymatgen.analysis.interfaces.coherent_interfaces import CoherentInterfaceBuilder, ZSLGenerator, Interface
from ase.build.tools import niggli_reduce

from ..convert import to_ase, from_ase, to_pymatgen, from_pymatgen
from .slab import BaseSlabConfiguration, SlabConfiguration
from ...material import Material

TerminationPair = Tuple[str, str]
InterfacesType = List[Interface]
InterfacesDataType = Dict[Tuple, List[Interface]]


class StrainModes(str, Enum):
    strain = "strain"
    von_mises_strain = "von_mises_strain"
    mean_abs_strain = "mean_abs_strain"


class ZSLStrainMatchingParameters(BaseModel):
    max_area: float = 50.0
    max_area_ratio_tol: float = 0.09
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
    ):
        super().__init__()
        self.substrate_configuration = substrate_configuration
        self.film_configuration = film_configuration
        self.termination_pair = termination_pair
        self.vacuum: float = vacuum
        self.distance_z: float = distance_z

    @property
    def bulk(self):
        # TODO: Return the bulk structure of the interface configuration, with no vacuum
        return self.get_material()

    @property
    def miller_indices(self):
        return (0, 0, 1)

    def get_material(self, scale_film_to_fit: bool = False):
        interface_builder = SimpleInterfaceBuilder(
            strain_matching_parameters=NoStrainMatchingParameters(scale_film_to_fit=scale_film_to_fit)
        )
        return interface_builder.get_material(self)


class NoStrainMatchingParameters(BaseModel):
    scale_film: bool = True


class BaseInterfaceBuilder:
    def __init__(
        self,
        strain_matching_parameters: Optional[Union[NoStrainMatchingParameters, ZSLStrainMatchingParameters]] = None,
    ):
        self.strain_matching_parameters = strain_matching_parameters
        pass

    def get_materials(self, interface_configuration: InterfaceConfiguration) -> List[Material]:
        raise NotImplementedError

    def get_material(self, interface_configuration: InterfaceConfiguration) -> Material:
        raise NotImplementedError


class SimpleInterfaceBuilder(BaseInterfaceBuilder):
    """
    Creates matching interface between substrate and film by straining the film to match the substrate.
    """

    def __init__(self, strain_matching_parameters: Optional[NoStrainMatchingParameters] = None):
        super().__init__(strain_matching_parameters=strain_matching_parameters)

    def get_material(self, interface_configuration: InterfaceConfiguration) -> Material:
        vacuum = interface_configuration.vacuum
        substrate_slab = interface_configuration.substrate_configuration.get_material(
            termination=interface_configuration.termination_pair[1]
        )
        film_slab = interface_configuration.film_configuration.get_material(
            termination=interface_configuration.termination_pair[0]
        )

        substrate_slab_ase = to_ase(substrate_slab)
        film_slab_ase = to_ase(film_slab)

        # Reduce the cell to Niggli form for correct alignment
        niggli_reduce(substrate_slab_ase)
        niggli_reduce(film_slab_ase)

        if self.strain_matching_parameters.scale_film:
            film_slab_ase.set_cell(substrate_slab_ase.cell, scale_atoms=True)
            film_slab_ase.wrap()

        # Calculate z-shift based on the new positions after scaling
        max_z_substrate = max(substrate_slab_ase.positions[:, 2])
        min_z_film = min(film_slab_ase.positions[:, 2])
        shift_z = max_z_substrate - min_z_film + interface_configuration.distance_z

        film_slab_ase.translate([0, 0, shift_z])

        interface_ase = substrate_slab_ase + film_slab_ase

        # Adjust the cell height to include the vacuum space above the film
        cell_c_with_vacuum = max(interface_ase.positions[:, 2]) + vacuum
        interface_ase.cell[2, 2] = cell_c_with_vacuum

        material_dict = from_ase(interface_ase)
        return Material(material_dict)


class ZSLInterfaceBuilder(BaseInterfaceBuilder):
    """
    Creates matching interface between substrate and film using the ZSL algorithm.
    """

    def __init__(self, strain_matching_parameters: ZSLStrainMatchingParameters):
        super().__init__(strain_matching_parameters=strain_matching_parameters)
        self.strain_matching_parameters = strain_matching_parameters.dict()

    def get_material(self, interface_configuration: InterfaceConfiguration) -> Material:
        interface = self.get_sorted_interfaces()[0]
        interface_ase = to_ase(interface)
        material_dict = from_ase(interface_ase)
        return Material(material_dict)

    def get_sorted_interfaces(self):
        interfaces = self.get_interface_structures()
        interfaces = sorted(interfaces, key=lambda x: np.mean(np.abs(x.interface_properties["strain"])))
        return [Material(from_pymatgen(interface)) for interface in interfaces]

    def get_interface_structures(self, interface_configuration: InterfaceConfiguration) -> InterfacesType:
        interfaces = self.generator(interface_configuration).get_interfaces(
            termination=interface_configuration.termination_pair,
            gap=interface_configuration.distance_z,
            film_thickness=interface_configuration.film_configuration.thickness,
            substrate_thickness=interface_configuration.substrate_configuration.thickness,
            in_layers=True,
        )
        return list(interfaces)

    def generator(self, interface_configuration: InterfaceConfiguration) -> CoherentInterfaceBuilder:
        generator = ZSLGenerator(**self.strain_matching_parameters)
        return CoherentInterfaceBuilder(
            substrate_structure=to_pymatgen(interface_configuration.substrate_configuration.bulk),
            film_structure=to_pymatgen(interface_configuration.film_configuration.bulk),
            substrate_miller=interface_configuration.substrate_configuration.miller_indices,
            film_miller=interface_configuration.film_configuration.miller_indices,
            zslgen=generator,
        )
