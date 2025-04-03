from typing import Any, List, Optional

import numpy as np
from ase.build.tools import niggli_reduce
from mat3ra.made.lattice import Lattice
from mat3ra.made.material import Material
from mat3ra.made.utils import create_2d_supercell_matrices, get_angle_from_rotation_matrix_2d
from pydantic import BaseModel
from pymatgen.analysis.interfaces.coherent_interfaces import (
    CoherentInterfaceBuilder,
    ZSLGenerator,
)

from .commensurate_lattice_pair import CommensurateLatticePair
from .configuration import (
    InterfaceConfiguration,
    NanoRibbonTwistedInterfaceConfiguration,
    TwistedInterfaceConfiguration,
)
from .enums import StrainModes, angle_to_supercell_matrix_values_for_hex
from .termination_pair import TerminationPair, safely_select_termination_pair
from .utils import interface_patch_with_mean_abs_strain, remove_duplicate_interfaces
from ..mixins import (
    ConvertGeneratedItemsASEAtomsMixin,
    ConvertGeneratedItemsPymatgenStructureMixin,
)
from ..nanoribbon import NanoribbonConfiguration, create_nanoribbon
from ..slab import create_slab, Termination, SlabConfiguration
from ..supercell import create_supercell
from ..utils import merge_materials
from ...analyze.other import get_chemical_formula
from ...build import BaseBuilder, BaseBuilderParameters
from ...convert import to_ase, from_ase, to_pymatgen, PymatgenInterface, ASEAtoms
from ...modify import (
    translate_to_z_level,
    rotate,
    translate_by_vector,
    add_vacuum_sides,
    add_vacuum,
)


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
    scale_film: bool = True  # Whether to scale the film to match the substrate
    create_slabs: bool = True  # Whether to create slabs from the configurations or use the bulk


class SimpleInterfaceBuilder(ConvertGeneratedItemsASEAtomsMixin, InterfaceBuilder):
    """
    Creates matching interface between substrate and film by straining the film to match the substrate.
    """

    _BuildParametersType = SimpleInterfaceBuilderParameters
    _DefaultBuildParameters = SimpleInterfaceBuilderParameters(scale_film=True)
    _GeneratedItemType: type(ASEAtoms) = ASEAtoms  # type: ignore

    def __preprocess_slab_configuration(
        self, configuration: SlabConfiguration, termination: Termination, create_slabs=False
    ) -> ASEAtoms:
        slab = create_slab(configuration, termination) if create_slabs else configuration.bulk
        ase_slab = to_ase(slab)

        niggli_reduce(ase_slab)
        return ase_slab

    @staticmethod
    def __combine_two_slabs_ase(substrate_slab_ase: ASEAtoms, film_slab_ase: ASEAtoms, distance_z: float) -> ASEAtoms:
        total_z_height = substrate_slab_ase.cell[2][2] + film_slab_ase.cell[2][2] + distance_z
        substrate_slab_ase.cell[2][2] = total_z_height
        film_slab_ase.cell[2][2] = total_z_height
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
            configuration.film_configuration,
            configuration.film_termination,
            create_slabs=self.build_parameters.create_slabs,
        )
        substrate_slab_ase = self.__preprocess_slab_configuration(
            configuration.substrate_configuration,
            configuration.substrate_termination,
            create_slabs=self.build_parameters.create_slabs,
        )

        if self.build_parameters.scale_film:
            temp_cell = [substrate_slab_ase.cell[0], substrate_slab_ase.cell[1], film_slab_ase.cell[2]]
            # We scale the film in XY direction only, and use two scaling steps and a temporary cell
            film_slab_ase.set_cell(temp_cell, scale_atoms=True)
            film_slab_ase.set_cell(substrate_slab_ase.cell, scale_atoms=False)
            film_slab_ase.wrap()

        interface_ase = self.__combine_two_slabs_ase(substrate_slab_ase, film_slab_ase, configuration.distance_z)
        interface_ase_with_vacuum = self.__add_vacuum_along_c_ase(interface_ase, configuration.vacuum)

        return [interface_ase_with_vacuum]

    def _post_process(self, items: List[_GeneratedItemType], post_process_parameters=None) -> List[Material]:
        return [Material.create(from_ase(slab)) for slab in items]


########################################################################################
#                       Strain Matching Interface Builders                             #
########################################################################################
class StrainMatchingInterfaceBuilderParameters(BaseBuilderParameters):
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
    strain_matching_parameters: ZSLStrainMatchingParameters = ZSLStrainMatchingParameters()


class ZSLStrainMatchingInterfaceBuilder(ConvertGeneratedItemsPymatgenStructureMixin, StrainMatchingInterfaceBuilder):
    """
    Creates matching interface between substrate and film using the ZSL algorithm.
    """

    _BuildParametersType: type(  # type: ignore
        ZSLStrainMatchingInterfaceBuilderParameters
    ) = ZSLStrainMatchingInterfaceBuilderParameters  # type: ignore
    _GeneratedItemType: PymatgenInterface = PymatgenInterface  # type: ignore

    def _generate(self, configuration: InterfaceConfiguration) -> List[PymatgenInterface]:
        generator = ZSLGenerator(**self.build_parameters.strain_matching_parameters.model_dump())
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


class NanoRibbonTwistedInterfaceBuilder(BaseBuilder):
    _GeneratedItemType = Material
    _ConfigurationType: type(  # type: ignore
        NanoRibbonTwistedInterfaceConfiguration
    ) = NanoRibbonTwistedInterfaceConfiguration  # type: ignore

    def _generate(self, configuration: _ConfigurationType) -> List[Material]:
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
        top_ribbon = rotate(top_ribbon, [0, 0, 1], configuration.twist_angle, wrap=False)

        translation_vector = [0, 0, configuration.distance_z]
        top_ribbon = translate_by_vector(top_ribbon, translation_vector, use_cartesian_coordinates=True)

        merged_material = merge_materials([bottom_ribbon, top_ribbon])
        merged_material_vacuum_x = add_vacuum_sides(merged_material, configuration.vacuum_x, on_x=True)
        merged_material_vacuum_xy = add_vacuum_sides(merged_material_vacuum_x, configuration.vacuum_y, on_y=True)

        return [merged_material_vacuum_xy]

    def _update_material_name(
        self, material: Material, configuration: NanoRibbonTwistedInterfaceConfiguration
    ) -> Material:
        material.name = f"Twisted Nanoribbon Interface ({configuration.twist_angle:.2f} degrees)"
        return material


class CommensurateLatticeTwistedInterfaceBuilderParameters(BaseModel):
    """
    Parameters for the commensurate lattice interface builder.

    Args:
        max_supercell_matrix_int (Optional[int]): The maximum integer for the transformation matrices.
            If not provided, it will be determined based on the target angle and the lattice vectors automatically.
        limit_max_int (Optional[int]): The limit for the maximum integer for the transformation matrices when searching
        angle_tolerance (float): The tolerance for the angle between the commensurate lattices
            and the target angle, in degrees.
        return_first_match (bool): Whether to return the first match or all matches.
    """

    max_supercell_matrix_int: Optional[int] = None
    limit_max_int: Optional[int] = 42
    angle_tolerance: float = 0.1
    return_first_match: bool = False


class CommensurateLatticeTwistedInterfaceBuilder(BaseBuilder):
    _GeneratedItemType: type(CommensurateLatticePair) = CommensurateLatticePair  # type: ignore
    _ConfigurationType: type(TwistedInterfaceConfiguration) = TwistedInterfaceConfiguration  # type: ignore

    def _generate(self, configuration: _ConfigurationType) -> List[_GeneratedItemType]:
        film = configuration.film
        # substrate = configuration.substrate
        a = film.lattice.vector_arrays[0][:2]
        b = film.lattice.vector_arrays[1][:2]
        max_int = self.build_parameters.max_supercell_matrix_int or self.__get_initial_guess_for_max_int(
            film, configuration.twist_angle
        )
        commensurate_lattice_pairs: List[CommensurateLatticePair] = []
        while not commensurate_lattice_pairs and max_int < self.build_parameters.limit_max_int:
            print(f"Trying max_int = {max_int}")
            commensurate_lattice_pairs = self.__generate_commensurate_lattices(
                configuration, a, b, max_int, configuration.twist_angle
            )
            max_int += 1

        return commensurate_lattice_pairs

    def __get_initial_guess_for_max_int(self, film, target_angle: float) -> int:
        """
        Determine the maximum integer for the transformation matrices based on the target angle.

        Args:
            a (List[float]): The a lattice vector.
            b (List[float]): The b lattice vector.
            target_angle (float): The target twist angle, in degrees.

        Returns:
            int: The maximum integer for the transformation matrices.
        """
        if film.lattice.type == Lattice.__types__.HEX:
            # getting max int of the matrix that has angle closest to target angle
            xy_supercell_matrix_for_closest_angle = min(
                angle_to_supercell_matrix_values_for_hex, key=lambda x: abs(x["angle"] - target_angle)
            )
            # Get maximum absolute value from the supercell matrix values
            return max(abs(x) for row in xy_supercell_matrix_for_closest_angle["xy_supercell"] for x in row)
        return 1

    def __generate_commensurate_lattices(
        self,
        configuration: TwistedInterfaceConfiguration,
        a: List[float],
        b: List[float],
        max_supercell_matrix_element_int: int,
        target_angle: float = 0.0,
    ) -> List[CommensurateLatticePair]:
        """
        Generate all commensurate lattices for a given target angle and filter by closeness to target angle.

        Args:
            configuration (TwistedInterfaceConfiguration): The configuration for the twisted interface.
            a (List[float]): The a lattice vector.
            b (List[float]): The b lattice vector.
            max_supercell_matrix_element_int (int): The maximum integer for the transformation matrices.
            target_angle (float): The target twist angle, in degrees.

        Returns:
            List[CommensurateLatticePair]: The list of commensurate lattice pairs
        """
        # Generate the supercell matrices using the calculated max_supercell_matrix_element_int
        matrices = create_2d_supercell_matrices(max_supercell_matrix_element_int)
        matrix_ab = np.array([a, b])
        matrix_ab_inverse = np.linalg.inv(matrix_ab)

        solutions: List[CommensurateLatticePair] = []
        for index1, matrix1 in enumerate(matrices):
            for index2, matrix2 in enumerate(matrices[0 : index1 + 1]):
                matrix2_inverse = np.linalg.inv(matrix2)
                intermediate_product = matrix2_inverse @ matrix1
                product = matrix_ab_inverse @ intermediate_product @ matrix_ab
                angle = get_angle_from_rotation_matrix_2d(product)
                if angle is not None:
                    size_metric = np.linalg.det(matrix_ab_inverse @ matrix1 @ matrix_ab)

                    if np.abs(angle - target_angle) < self.build_parameters.angle_tolerance:
                        print(f"Found commensurate lattice with angle {angle} and size metric {size_metric}")
                        solutions.append(
                            CommensurateLatticePair(
                                configuration=configuration,
                                matrix1=matrix1,
                                matrix2=matrix2,
                                angle=angle,
                                size_metric=size_metric,
                            )
                        )
                        if self.build_parameters.return_first_match:
                            return solutions
                else:
                    continue
        return solutions

    def _post_process(
        self,
        items: List[_GeneratedItemType],
        post_process_parameters=None,
    ) -> List[Material]:
        interfaces = []
        for item in items:
            new_substrate = create_supercell(item.configuration.film, item.matrix1.tolist())
            new_film = create_supercell(item.configuration.substrate, item.matrix2.tolist())
            new_film = translate_by_vector(
                new_film, [0, 0, item.configuration.distance_z], use_cartesian_coordinates=True
            )
            interface = merge_materials([new_substrate, new_film], merge_dangerously=True)
            interface.metadata["actual_twist_angle"] = item.angle
            if item.configuration.vacuum != 0:
                interface = add_vacuum(interface, item.configuration.vacuum)
            interfaces.append(interface)
        return interfaces

    def _update_material_metadata(self, material, configuration) -> Material:
        updated_material = super()._update_material_metadata(material, configuration)
        if "actual_twist_angle" in material.metadata:
            updated_material.metadata["build"]["configuration"]["actual_twist_angle"] = material.metadata[
                "actual_twist_angle"
            ]
        return updated_material

    def _update_material_name(
        self, material: Material, configuration: NanoRibbonTwistedInterfaceConfiguration
    ) -> Material:
        material.name = f"Twisted Bilayer Interface ({material.metadata['actual_twist_angle']:.2f} degrees)"
        return material
