from typing import Any, List, Optional

import numpy as np
from mat3ra.made.tools.build.supercell import create_supercell

from ..utils import merge_materials
from ...modify import (
    translate_to_z_level,
    rotate_material,
    translate_by_vector,
    add_vacuum_sides,
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
    Configuration for creating a twisted interface between two nano ribbons with specified twist angle.

    Args:
        film (Material): The film material.
        substrate (Material): The substrate material.
        twist_angle (float): Twist angle in degrees.
        ribbon_width (int): Width of the nanoribbon in unit cells.
        ribbon_length (int): Length of the nanoribbon in unit cells.
        distance_z (float): Vertical distance between layers in Angstroms.
        vacuum_x (float): Vacuum along x on both sides, in Angstroms.
        vacuum_y (float): Vacuum along y on both sides, in Angstroms.
    """

    ribbon_width: int = 1
    ribbon_length: int = 1
    vacuum_x: float = 5.0
    vacuum_y: float = 5.0

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
        )
        bottom_ribbon = create_nanoribbon(bottom_nanoribbon_configuration)
        top_ribbon_configuration = NanoribbonConfiguration(
            material=configuration.film,
            width=configuration.ribbon_width,
            length=configuration.ribbon_length,
        )
        top_ribbon = create_nanoribbon(top_ribbon_configuration)
        top_ribbon = rotate_material(top_ribbon, [0, 0, 1], configuration.twist_angle, wrap=False)

        translation_vector = [0, 0, configuration.distance_z]
        top_ribbon = translate_by_vector(top_ribbon, translation_vector, use_cartesian_coordinates=True)

        merged_material = merge_materials([bottom_ribbon, top_ribbon])
        merged_material_vacuum_x = add_vacuum_sides(merged_material, configuration.vacuum_x, on_x=True)
        merged_material_vacuum_xy = add_vacuum_sides(merged_material_vacuum_x, configuration.vacuum_y, on_y=True)

        return [merged_material_vacuum_xy]

    def _update_material_name(
        self, material: Material, configuration: NanoRibbonTwistedInterfaceConfiguration
    ) -> Material:
        material.name = f"Twisted Nanoribbon Interface ({configuration.twist_angle:.2f}°)"
        return material


class CommensurateLatticeInterfaceBuilderParameters(BaseModel):
    max_search: int = 10
    angle_tolerance: float = 1.0


class CommensurateLatticeInterfaceBuilder(BaseBuilder):
    _GeneratedItemType = Material
    _ConfigurationType = TwistedInterfaceConfiguration

    def _generate(self, configuration: TwistedInterfaceConfiguration) -> List[Material]:
        film = configuration.film
        substrate = configuration.substrate
        max_search = self.build_parameters.max_search
        a1 = film.lattice.vector_arrays[0][:2]
        a2 = film.lattice.vector_arrays[1][:2]
        commensurate_lattices = self.__generate_commensurate_lattices(a1, a2, max_search)
        commensurate_lattices_with_angle = [
            lattice
            for lattice in commensurate_lattices
            if np.abs(lattice["angle"] - configuration.twist_angle) < self.build_parameters.angle_tolerance
        ]
        interfaces = []
        for lattice in commensurate_lattices_with_angle:
            new_substrate = create_supercell(film, lattice["matrix1"].tolist())
            new_film = create_supercell(substrate, lattice["matrix2"].tolist())
            new_film = translate_by_vector(new_film, [0, 0, configuration.distance_z], use_cartesian_coordinates=True)
            try:
                interface = merge_materials([new_substrate, new_film])
                interfaces.append(interface)
            except Exception as e:
                print(e)
        return interfaces

    @staticmethod
    def __create_matrices(max_search: int):
        matrices = []
        for s11 in range(-max_search, max_search + 1):
            for s12 in range(-max_search, max_search + 1):
                for s21 in range(-max_search, max_search + 1):
                    for s22 in range(-max_search, max_search + 1):
                        # Non-zero area constraint
                        matrix = np.array([[s11, s12], [s21, s22]])
                        determinant = np.linalg.det(matrix)
                        if determinant == 0 or determinant < 0:
                            continue
                        matrices.append(matrix)
        return matrices

    @staticmethod
    def __solve_angle_from_rotation_matrix(matrix, zero_tolerance=1e-6, round_digits=3):
        if matrix.shape != (2, 2):
            raise ValueError("Input matrix must be 2x2")
        if np.abs(np.linalg.det(matrix) - 1) > zero_tolerance:
            raise ValueError("Matrix must be orthogonal (determinant = 1)")
        if not np.all(np.abs(matrix) <= 1):
            raise ValueError("Matrix have all elements less than 1")
        cos_theta = matrix[0, 0]
        sin_theta = matrix[1, 0]
        angle_rad = np.arctan2(sin_theta, cos_theta)
        angle_deg = np.round(np.degrees(angle_rad), round_digits)

        return angle_deg

    def __generate_commensurate_lattices(self, a1: List[float], a2: List[float], max_search: int = 10):
        """
        Generate all commensurate lattices for a given search range.
        """
        matrices = self.__create_matrices(max_search)
        matrix_a1a2 = np.array([a1, a2])
        matrix_a1a2_inverse = np.linalg.inv(matrix_a1a2)

        solutions = []
        for index1, matrix1 in enumerate(matrices):
            for index2, matrix2 in enumerate(matrices[0 : index1 + 1]):
                matrix2_inverse = np.linalg.inv(matrix2)
                intermediate_product = matrix2_inverse @ matrix1
                product = matrix_a1a2_inverse @ intermediate_product @ matrix_a1a2
                try:
                    angle = self.__solve_angle_from_rotation_matrix(product)
                    size_metric = np.linalg.det(matrix_a1a2_inverse @ matrix1 @ matrix_a1a2)
                    solutions.append(
                        {"matrix1": matrix1, "matrix2": matrix2, "angle": angle, "size_metric": size_metric}
                    )
                except ValueError:
                    continue
        return solutions
