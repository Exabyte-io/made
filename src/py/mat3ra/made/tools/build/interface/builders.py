from typing import Any, List, Optional

import numpy as np
from ..supercell import create_supercell
from ..utils import merge_materials
from ...modify import (
    translate_to_z_level,
    rotate_material,
    translate_by_vector,
)
from pydantic import BaseModel, Field
from ase.build.tools import niggli_reduce
from pymatgen.analysis.interfaces.coherent_interfaces import (
    CoherentInterfaceBuilder,
    ZSLGenerator,
)

from ..nanoribbon import NanoribbonConfiguration, create_nanoribbon

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
from ...convert import to_ase, from_ase, to_pymatgen, PymatgenInterface, ASEAtoms
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
    twist_angle: float = Field(0, description="Twist angle in degrees")
    distance_z: float = 3.0

    @property
    def _json(self):
        return {
            "type": self.get_cls_name(),
            "film": self.film.to_json(),
            "substrate": self.substrate.to_json(),
            "twist_angle": self.twist_angle,
            "distance_z": self.distance_z,
        }


class NanoRibbonTwistedInterfaceConfiguration(TwistedInterfaceConfiguration):
    """
    Configuration for creating a twisted interface between two materials using nanoribbons.

    Args:
        film (Material): The film material.
        substrate (Material): The substrate material.
        twist_angle (float): Twist angle in degrees.
        ribbon_width (int): Width of the nanoribbon in unit cells.
        ribbon_length (int): Length of the nanoribbon in unit cells.
        distance_z (float): Vertical distance between layers in Angstroms.
        vacuum_x (int): Vacuum padding in x direction in unit cells.
        vacuum_y (int): Vacuum padding in y direction in unit cells.
    """

    ribbon_width: int = 1
    ribbon_length: int = 1
    vacuum_x: int = 2
    vacuum_y: int = 2

    @property
    def _json(self):
        return {
            **super()._json,
            "type": self.get_cls_name(),
            "ribbon_width": self.ribbon_width,
            "ribbon_length": self.ribbon_length,
            "vacuum_x": self.vacuum_x,
            "vacuum_y": self.vacuum_y,
        }


class NanoRibbonTwistedInterfaceBuilder(BaseBuilder):
    _GeneratedItemType = Material
    _ConfigurationType = NanoRibbonTwistedInterfaceConfiguration

    def _generate(self, configuration: NanoRibbonTwistedInterfaceConfiguration) -> List[Material]:
        bottom_nanoribbon_configuration = NanoribbonConfiguration(
            material=configuration.substrate,
            width=configuration.ribbon_width,
            length=configuration.ribbon_length,
            vacuum_width=configuration.vacuum_x,
            vacuum_length=configuration.vacuum_y,
        )
        bottom_ribbon = create_nanoribbon(bottom_nanoribbon_configuration)

        top_ribbon_configuration = NanoribbonConfiguration(
            material=configuration.film,
            width=configuration.ribbon_width,
            length=configuration.ribbon_length,
            vacuum_width=configuration.vacuum_x,
            vacuum_length=configuration.vacuum_y,
        )
        top_ribbon = create_nanoribbon(top_ribbon_configuration)
        top_ribbon = rotate_material(top_ribbon, [0, 0, 1], configuration.twist_angle, wrap=False)

        translation_vector = [0, 0, configuration.distance_z]
        top_ribbon = translate_by_vector(top_ribbon, translation_vector, use_cartesian_coordinates=True)
        merged_material = merge_materials([bottom_ribbon, top_ribbon])

        return [merged_material]

    def _update_material_name(
        self, material: Material, configuration: NanoRibbonTwistedInterfaceConfiguration
    ) -> Material:
        material.name = f"Twisted Nanoribbon Interface ({configuration.twist_angle:.2f}°)"
        return material


class CommensurateSuperCellTwistedInterfaceConfiguration(TwistedInterfaceConfiguration):
    max_supercell_size: int = Field(10, description="Maximum supercell size to consider")

    @property
    def _json(self):
        return {
            **super()._json,
            "type": self.get_cls_name(),
            "max_supercell_size": self.max_supercell_size,
        }


class CommensurateSuperCellTwistedInterfaceBuilder(BaseBuilder):
    _GeneratedItemType = Material
    _ConfigurationType = CommensurateSuperCellTwistedInterfaceConfiguration

    def _generate(self, configuration: CommensurateSuperCellTwistedInterfaceConfiguration) -> List[Material]:
        M1, M2, strain = self.find_csl_matrices(
            configuration.twist_angle,
            max_index=configuration.max_supercell_size,
            angle_tolerance=0.5,
            strain_tolerance=0.05,
        )

        film_1 = create_supercell(configuration.substrate, M1)
        film_2 = create_supercell(configuration.film, M2)
        film_2 = translate_by_vector(film_2, [0, 0, configuration.distance_z], use_cartesian_coordinates=True)

        bilayer = merge_materials([film_2, film_1], merge_dangerously=True)
        return [bilayer]

    @staticmethod
    def __rotation_matrix_2d(angle: float) -> np.ndarray:
        """Create a 2D rotation matrix for the given angle in degrees."""
        theta = np.radians(angle)
        return np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])

    def find_csl_matrices(
        self, angle: float, max_index: int = 10, angle_tolerance: float = 0.1, strain_tolerance: float = 0.01
    ):
        R_target = self.__rotation_matrix_2d(angle)
        best_strain = float("inf")
        best_M1 = best_M2 = None

        for i in range(1, max_index + 1):
            for j in range(max_index + 1):
                M1 = np.array([[i, -j], [j, i]])
                M2_approx = R_target @ M1
                M2 = np.round(M2_approx).astype(int)
                R_actual = M2 @ np.linalg.inv(M1)

                angle_actual = np.degrees(np.arctan2(R_actual[1, 0], R_actual[0, 0]))
                if abs(angle_actual - angle) > angle_tolerance:
                    continue

                strain = np.linalg.norm(M2_approx - M2) / np.linalg.norm(M2_approx)
                if strain < strain_tolerance and strain < best_strain:
                    best_strain = strain
                    best_M1 = np.array([[i, -j, 0], [j, i, 0], [0, 0, 1]])
                    best_M2 = np.array([M2[0, 0], M2[0, 1], 0, M2[1, 0], M2[1, 1], 0, 0, 0, 1]).reshape(3, 3)

        if best_M1 is None:
            raise ValueError(f"No approximate CSL found for angle {angle} within tolerances")

        return best_M1, best_M2, best_strain

    def _update_material_name(
        self, material: Material, configuration: CommensurateSuperCellTwistedInterfaceConfiguration
    ) -> Material:
        material.name = f"Commensurate Twisted Interface : {material.metadata['build']['actual_angle']:.2f} degrees"
        return material
