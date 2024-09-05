from typing import Any, List, Optional, Tuple

import numpy as np
from mat3ra.made.tools.build.supercell import create_supercell
from mat3ra.made.tools.build.utils import merge_materials
from mat3ra.made.tools.modify import translate_to_z_level, rotate_material, translate_by_vector, filter_by_box
from pydantic import BaseModel, Field
from ase.build.tools import niggli_reduce
from pymatgen.analysis.interfaces.coherent_interfaces import (
    CoherentInterfaceBuilder,
    ZSLGenerator,
)
from ...third_party import PymatgenStructure

from mat3ra.made.material import Material
from .enums import StrainModes
from .configuration import InterfaceConfiguration
from .termination_pair import TerminationPair, safely_select_termination_pair
from .utils import interface_patch_with_mean_abs_strain, remove_duplicate_interfaces
from ..mixins import (
    ConvertGeneratedItemsASEAtomsMixin,
    ConvertGeneratedItemsPymatgenStructureMixin,
)
from ..slab import create_slab, Termination
from ..slab.configuration import SlabConfiguration
from ...analyze import get_chemical_formula
from ...convert import to_ase, from_ase, to_pymatgen, PymatgenInterface, ASEAtoms, from_pymatgen
from ...build import BaseBuilder, BaseConfiguration


class InterfaceBuilderParameters(BaseModel):
    pass


class InterfaceBuilder(BaseBuilder):
    _BuildParametersType = InterfaceBuilderParameters
    _ConfigurationType: type(InterfaceConfiguration) = InterfaceConfiguration  # type: ignore

    def _update_material_name(self, material: Material, configuration: InterfaceConfiguration) -> Material:
        film_formula = get_chemical_formula(configuration.film_configuration.bulk)
        substrate_formula = get_chemical_formula(configuration.substrate_configuration.bulk)
        film_miller_indices = "".join([str(i) for i in configuration.film_configuration.miller_indices])
        substrate_miller_indices = "".join([str(i) for i in configuration.substrate_configuration.miller_indices])
        new_name = f"{film_formula}({film_miller_indices})-{substrate_formula}({substrate_miller_indices}), Interface"
        material.name = new_name
        return material


########################################################################################
#                           Simple Interface Builder                                   #
########################################################################################
class SimpleInterfaceBuilderParameters(InterfaceBuilderParameters):
    scale_film: bool = True


class SimpleInterfaceBuilder(ConvertGeneratedItemsASEAtomsMixin, InterfaceBuilder):
    """
    Creates matching interface between substrate and film by straining the film to match the substrate.
    """

    _BuildParametersType = SimpleInterfaceBuilderParameters
    _DefaultBuildParameters = SimpleInterfaceBuilderParameters(scale_film=True)
    _GeneratedItemType: type(ASEAtoms) = ASEAtoms  # type: ignore

    @staticmethod
    def __preprocess_slab_configuration(configuration: SlabConfiguration, termination: Termination):
        slab = create_slab(configuration, termination)
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


########################################################################################
#                       Strain Matching Interface Builders                             #
########################################################################################
class StrainMatchingInterfaceBuilderParameters(BaseModel):
    strain_matching_parameters: Optional[Any] = None


class StrainMatchingInterfaceBuilder(InterfaceBuilder):
    _BuildParametersType = StrainMatchingInterfaceBuilderParameters  # type: ignore

    def _update_material_name(self, material: Material, configuration: InterfaceConfiguration) -> Material:
        updated_material = super()._update_material_name(material, configuration)
        if StrainModes.mean_abs_strain in material.metadata:
            strain = material.metadata[StrainModes.mean_abs_strain]
            new_name = f"{updated_material.name}, Strain {strain*100:.3f}pct"
            updated_material.name = new_name
        return updated_material


class ZSLStrainMatchingParameters(BaseModel):
    max_area: float = 50.0
    max_area_ratio_tol: float = 0.09
    max_length_tol: float = 0.03
    max_angle_tol: float = 0.01


class ZSLStrainMatchingInterfaceBuilderParameters(StrainMatchingInterfaceBuilderParameters):
    strain_matching_parameters: ZSLStrainMatchingParameters


class ZSLStrainMatchingInterfaceBuilder(ConvertGeneratedItemsPymatgenStructureMixin, StrainMatchingInterfaceBuilder):
    """
    Creates matching interface between substrate and film using the ZSL algorithm.
    """

    _BuildParametersType: type(  # type: ignore
        ZSLStrainMatchingInterfaceBuilderParameters
    ) = ZSLStrainMatchingInterfaceBuilderParameters  # type: ignore
    _GeneratedItemType: PymatgenInterface = PymatgenInterface  # type: ignore

    def _generate(self, configuration: InterfaceConfiguration) -> List[PymatgenInterface]:
        generator = ZSLGenerator(**self.build_parameters.strain_matching_parameters.dict())
        substrate_with_atoms_translated_to_bottom = translate_to_z_level(
            configuration.substrate_configuration.bulk, "bottom"
        )
        film_with_atoms_translated_to_bottom = translate_to_z_level(configuration.film_configuration.bulk, "bottom")
        builder = CoherentInterfaceBuilder(
            substrate_structure=to_pymatgen(substrate_with_atoms_translated_to_bottom),
            film_structure=to_pymatgen(film_with_atoms_translated_to_bottom),
            substrate_miller=configuration.substrate_configuration.miller_indices,
            film_miller=configuration.film_configuration.miller_indices,
            zslgen=generator,
        )

        generated_termination_pairs = [
            TerminationPair.from_pymatgen(pymatgen_termination) for pymatgen_termination in builder.terminations
        ]
        termination_pair = safely_select_termination_pair(configuration.termination_pair, generated_termination_pairs)
        interfaces = builder.get_interfaces(
            termination=termination_pair.to_pymatgen(),
            gap=configuration.distance_z,
            vacuum_over_film=configuration.vacuum,
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


########################################################################################
#                       Twisted Interface Builders                                     #
########################################################################################


class TwistedInterfaceConfiguration(BaseConfiguration):
    film: Material
    substrate: Material
    twist_angle: float = Field(..., description="Twist angle in degrees")
    distance_z: float = 3.0

    @property
    def _json(self):
        json_data = super()._json
        json_data.update({"twist_angle": self.twist_angle})
        return json_data


class TwistedInterfaceBuilder(InterfaceBuilder):
    _ConfigurationType = TwistedInterfaceConfiguration

    def _generate(self, configuration: TwistedInterfaceConfiguration) -> List[Material]:
        film = configuration.film
        angle = configuration.twist_angle
        n = np.round(1 / np.sin(np.deg2rad(angle)))
        m = 2 * n
        distance_z_vector_crystal = film.basis.cell.convert_point_to_crystal([0, 0, configuration.distance_z])
        rotated_film = create_supercell(film, [[1, -1, 0], [1, 1, 0], [0, 0, 1]])
        rotated_film_scaled = create_supercell(rotated_film, scaling_factor=[n, m, 1])
        rotated_basis = rotate_material(rotated_film_scaled, [0, 0, 1], angle, wrap=False)
        rotated_basis_translated = translate_by_vector(rotated_basis, distance_z_vector_crystal)

        merged_material = merge_materials([rotated_film_scaled, rotated_basis_translated])

        merged_material_unrotated = rotate_material(merged_material, [0, 0, 1], -angle, wrap=False)
        merged_material_unrotated = filter_by_box(merged_material_unrotated, [0, 0, 0], [1, 1, 1])

        return [merged_material_unrotated]

    def _update_material_name(self, material: Material, configuration: TwistedInterfaceConfiguration) -> Material:
        film_formula = get_chemical_formula(configuration.film_configuration.bulk)
        substrate_formula = get_chemical_formula(configuration.substrate_configuration.bulk)
        film_miller = "".join(map(str, configuration.film_configuration.miller_indices))
        substrate_miller = "".join(map(str, configuration.substrate_configuration.miller_indices))

        new_name = f"{film_formula}({film_miller})-{substrate_formula}({substrate_miller}), Twist {configuration.twist_angle:.2f}°"
        if "mean_abs_strain" in material.metadata:
            strain = material.metadata["mean_abs_strain"]
            new_name += f", Strain {strain*100:.3f}%"

        material.name = new_name
        return material

    def _update_material_metadata(self, material: Material, configuration: TwistedInterfaceConfiguration) -> Material:
        material = super()._update_material_metadata(material, configuration)
        material.metadata["build"]["twist_angle"] = configuration.twist_angle
        return material

    def _post_process(self, items: List[PymatgenStructure], post_process_parameters=None) -> List[Material]:
        return [Material(from_pymatgen(structure)) for structure in items]
