from typing import Any, List, Optional

import numpy as np
from pydantic import BaseModel, Field
from ase.build.tools import niggli_reduce
from pymatgen.analysis.interfaces.coherent_interfaces import (
    CoherentInterfaceBuilder,
    ZSLGenerator,
)

from mat3ra.made.material import Material
from ....utils import create_2d_supercell_matrices, get_angle_from_rotation_matrix
from ...modify import (
    translate_to_z_level,
    rotate_material,
    translate_by_vector,
    add_vacuum_sides,
)
from ...analyze import get_chemical_formula
from ...convert import to_ase, from_ase, to_pymatgen, PymatgenInterface, ASEAtoms
from ...build import BaseBuilder, BaseConfiguration
from ..nanoribbon import NanoribbonConfiguration, create_nanoribbon
from ..supercell import create_supercell
from ..slab import create_slab, Termination, SlabConfiguration
from ..utils import merge_materials
from ..mixins import (
    ConvertGeneratedItemsASEAtomsMixin,
    ConvertGeneratedItemsPymatgenStructureMixin,
)

from .enums import StrainModes
from .configuration import InterfaceConfiguration
from .termination_pair import TerminationPair, safely_select_termination_pair
from .utils import interface_patch_with_mean_abs_strain, remove_duplicate_interfaces


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


class CommensurateLatticePair(BaseModel):
    class Config:
        arbitrary_types_allowed = True

    configuration: TwistedInterfaceConfiguration
    matrix1: np.ndarray
    matrix2: np.ndarray
    angle: float
    size_metric: float


class CommensurateLatticeInterfaceBuilderParameters(BaseModel):
    max_search: int = 10
    angle_tolerance: float = 0.1
    return_first_match: bool = False


class CommensurateLatticeInterfaceBuilder(BaseBuilder):
    _GeneratedItemType: type(CommensurateLatticePair) = CommensurateLatticePair  # type: ignore
    _ConfigurationType: type(TwistedInterfaceConfiguration) = TwistedInterfaceConfiguration  # type: ignore

    def _generate(self, configuration: _ConfigurationType) -> List[_GeneratedItemType]:
        film = configuration.film
        # substrate = configuration.substrate
        max_search = self.build_parameters.max_search
        a = film.lattice.vector_arrays[0][:2]
        b = film.lattice.vector_arrays[1][:2]
        commensurate_lattices = self.__generate_commensurate_lattices(a, b, max_search, configuration.twist_angle)
        commensurate_lattice_pairs = [
            CommensurateLatticePair(configuration=configuration, **lattice) for lattice in commensurate_lattices
        ]
        return commensurate_lattice_pairs

    def __generate_commensurate_lattices(
        self, a1: List[float], a2: List[float], max_search: int = 10, target_angle: float = 0.0
    ):
        """
        Generate all commensurate lattices for a given search range.
        """
        matrices = create_2d_supercell_matrices(max_search)
        matrix_a1a2 = np.array([a1, a2])
        matrix_a1a2_inverse = np.linalg.inv(matrix_a1a2)

        solutions = []
        for index1, matrix1 in enumerate(matrices):
            for index2, matrix2 in enumerate(matrices[0 : index1 + 1]):
                matrix2_inverse = np.linalg.inv(matrix2)
                intermediate_product = matrix2_inverse @ matrix1
                product = matrix_a1a2_inverse @ intermediate_product @ matrix_a1a2
                try:
                    angle = get_angle_from_rotation_matrix(product)
                    size_metric = np.linalg.det(matrix_a1a2_inverse @ matrix1 @ matrix_a1a2)

                    if np.abs(angle - target_angle) < self.build_parameters.angle_tolerance:
                        print(f"Found commensurate lattice with angle {angle} and size metric {size_metric}")
                        solutions.append(
                            {"matrix1": matrix1, "matrix2": matrix2, "angle": angle, "size_metric": size_metric}
                        )
                        if self.build_parameters.return_first_match:
                            return solutions
                except ValueError:
                    continue
        return solutions

    def _post_process(
        self,
        items: List[_GeneratedItemType],
        post_process_parameters: Optional[BaseBuilder._PostProcessParametersType] = None,
    ) -> List[Material]:
        interfaces = []
        for item in items:
            new_substrate = create_supercell(item.configuration.film, item.matrix1.tolist())
            new_film = create_supercell(item.configuration.substrate, item.matrix2.tolist())
            new_film = translate_by_vector(
                new_film, [0, 0, item.configuration.distance_z], use_cartesian_coordinates=True
            )
            try:
                interface = merge_materials([new_substrate, new_film])
                interfaces.append(interface)
            except Exception as e:
                print(e)
        return interfaces
