from typing import Any, List, Optional

import numpy as np
from pydantic import BaseModel
from ase.build.tools import niggli_reduce
from pymatgen.analysis.interfaces.coherent_interfaces import (
    CoherentInterfaceBuilder,
    ZSLGenerator,
)

from mat3ra.made.material import Material
from .enums import StrainModes
from .configuration import InterfaceConfiguration
from .termination_pair import TerminationPair
from .utils import interface_patch_with_mean_abs_strain, remove_duplicate_interfaces
from ..mixins import (
    ConvertGeneratedItemsASEAtomsMixin,
    ConvertGeneratedItemsPymatgenStructureMixin,
)
from ..slab.configuration import SlabConfiguration
from ...convert import to_ase, from_ase, to_pymatgen, PymatgenInterface, ASEAtoms
from ...build import BaseBuilder


class InterfaceBuilderParameters(BaseModel):
    pass


class InterfaceBuilder(BaseBuilder):
    _BuildParametersType = InterfaceBuilderParameters
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


class SimpleInterfaceBuilderParameters(InterfaceBuilderParameters):
    scale_film: bool = True


class SimpleInterfaceBuilder(InterfaceBuilder, ConvertGeneratedItemsASEAtomsMixin):
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

    def _post_process(self, items: List[_GeneratedItemType], post_process_parameters=None) -> List[Material]:
        return [Material(from_ase(slab)) for slab in items]


class StrainMatchingInterfaceBuilderParameters(BaseModel):
    strain_matching_parameters: Optional[Any] = None


class StrainMatchingInterfaceBuilder(InterfaceBuilder):
    _BuildParametersType = StrainMatchingInterfaceBuilderParameters

    def _update_material_name(self, material: Material, configuration: InterfaceConfiguration) -> Material:
        updated_material = super()._update_material_name(material, configuration)
        if StrainModes.mean_abs_strain in material.metadata:
            strain = material.metadata[StrainModes.mean_abs_strain]
            new_name = f"{updated_material.name}, Strain (mean_abs_strain={strain:.3f}%)"
            updated_material.name = new_name
        return material


class ZSLStrainMatchingParameters(BaseModel):
    max_area: float = 50.0
    max_area_ratio_tol: float = 0.09
    max_length_tol: float = 0.03
    max_angle_tol: float = 0.01


class ZSLStrainMatchingInterfaceBuilderParameters(StrainMatchingInterfaceBuilderParameters):
    strain_matching_parameters: ZSLStrainMatchingParameters


class ZSLStrainMatchingInterfaceBuilder(StrainMatchingInterfaceBuilder, ConvertGeneratedItemsPymatgenStructureMixin):
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

        # TODO: REMOVE when found the correct solution
        # Workaround: compare provided and generated terminations to find complete match,
        # if full match between terminations isn't found, get terminations with first part the same, before `_`

        generated_terminations = [
            TerminationPair(pymatgen_termination) for pymatgen_termination in builder.terminations
        ]
        provided_termination_pair = configuration.termination_pair
        hotfix_termination_pair = configuration.termination_pair

        if provided_termination_pair not in generated_terminations:
            for termination_pair in generated_terminations:
                if (
                    termination_pair.film_termination.split("_")[0]
                    == provided_termination_pair.film_termination.split("_")[0]
                    and termination_pair.substrate_termination.split("_")[0]
                    == provided_termination_pair.substrate_termination.split("_")[0]
                ):
                    hotfix_termination_pair = termination_pair
                    print("Interface will be built with terminations: ", hotfix_termination_pair)

        interfaces = builder.get_interfaces(
            termination=hotfix_termination_pair.self,
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
        materials = super()._post_process(items, post_process_parameters)
        strains = [interface.interface_properties[StrainModes.mean_abs_strain] for interface in items]

        for material, strain in zip(materials, strains):
            material.metadata["mean_abs_strain"] = strain
        return materials
