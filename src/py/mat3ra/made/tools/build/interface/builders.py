from typing import Any, List, Optional, Tuple

import numpy as np
from ..supercell import create_supercell
from ..utils import merge_materials
from ...modify import (
    translate_to_z_level,
    rotate_material,
    translate_by_vector,
    filter_by_triangle_projection,
    filter_by_box,
)
from pydantic import BaseModel, Field
from ase.build.tools import niggli_reduce
from pymatgen.analysis.interfaces.coherent_interfaces import (
    CoherentInterfaceBuilder,
    ZSLGenerator,
)

from ..nanoribbon import NanoribbonConfiguration, create_nanoribbon
from ...third_party import PymatgenStructure

from mat3ra.made.material import Material
from .enums import StrainModes, SupercellTypes
from .configuration import InterfaceConfiguration
from .termination_pair import TerminationPair, safely_select_termination_pair
from .utils import interface_patch_with_mean_abs_strain, remove_duplicate_interfaces
from ..mixins import (
    ConvertGeneratedItemsASEAtomsMixin,
    ConvertGeneratedItemsPymatgenStructureMixin,
)
from ..slab import create_slab, Termination
from ..slab.configuration import SlabConfiguration
from ...analyze import get_chemical_formula, calculate_moire_periodicity
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
    twist_angle: float = Field(0, description="Twist angle in degrees")
    distance_z: float = 3.0
    supercell_type: SupercellTypes = SupercellTypes.orthogonal

    @property
    def _json(self):
        json_data = super()._json
        json_data.update({"twist_angle": self.twist_angle})
        json_data.update({"distance_z": self.distance_z})
        json_data.update({"supercell_type": self.supercell_type})
        return json_data


class TwistedInterfaceBuilder(InterfaceBuilder):
    _ConfigurationType = TwistedInterfaceConfiguration

    def _generate(self, configuration: TwistedInterfaceConfiguration) -> List[Material]:
        film = configuration.film
        angle = configuration.twist_angle
        n, m = self._calculate_supercell_size(angle, 30)

        distance_z_vector_crystal = film.basis.cell.convert_point_to_crystal([0, 0, configuration.distance_z])
        rotated_film = create_supercell(film, [[1, -1, 0], [1, 1, 0], [0, 0, 1]])
        rotated_film_scaled = create_supercell(rotated_film, scaling_factor=[n, m, 1])
        rotated_basis = rotate_material(rotated_film_scaled, [0, 0, 1], angle, wrap=False)
        rotated_basis_translated = translate_by_vector(rotated_basis, distance_z_vector_crystal)

        merged_material = merge_materials([rotated_film_scaled, rotated_basis_translated])

        merged_material_unrotated = rotate_material(merged_material, [0, 0, 1], -angle, wrap=False)
        if configuration.supercell_type == SupercellTypes.orthogonal:
            final_material = self._create_orthogonal_supercell(rotated_film, angle, merged_material_unrotated)
        elif configuration.supercell_type == SupercellTypes.hexagonal:
            final_material = self._create_hex_supercell(film, angle, merged_material_unrotated)
        else:
            raise ValueError("Invalid supercell type")
        return [final_material]

    @staticmethod
    def _create_hex_supercell(film: Material, angle: float, merged_material: Material) -> Material:
        """
        Creates a hexagonal supercell for the twisted interface by filtering the atoms in the merged material
        in the shape of two triangles that form a parallelogram.

         Args:
            film (Material): The film material.
            angle (float): The twist angle in degrees.
            merged_material (Material): The merged material.

        Returns:
            Material: The final material with the hexagonal supercell.
        """
        p, q = calculate_moire_periodicity(film.lattice.a, film.lattice.b, angle)

        vertex_1 = [0.0, 0, 0]
        vertex_2 = [p / 2, np.sqrt(3) / 2 * q, 0]
        vertex_3 = [p, 0, 0]
        vertex_4 = [3 / 2 * p, np.sqrt(3) / 2 * q, 0]
        height = merged_material.lattice.c

        new_atoms_1 = filter_by_triangle_projection(
            merged_material,
            coordinate_1=vertex_1,
            coordinate_2=vertex_2,
            coordinate_3=vertex_3,
            max_z=height,
            use_cartesian_coordinates=True,
        )
        new_atoms_2 = filter_by_triangle_projection(
            merged_material,
            coordinate_1=vertex_2,
            coordinate_2=vertex_3,
            coordinate_3=vertex_4,
            max_z=height,
            use_cartesian_coordinates=True,
        )

        new_atoms = merge_materials([new_atoms_1, new_atoms_2])

        final_material = new_atoms.clone()
        final_material.set_new_lattice_vectors(vertex_2, vertex_3, [0, 0, height])

        return final_material

    @staticmethod
    def _create_orthogonal_supercell(film: Material, angle: float, merged_material: Material) -> Material:
        """
        Creates an orthogonal supercell for the twisted interface by filtering the atoms in the merged material
        in the shape of a rectangular box with the sides equal to the Moire pattern periodicity vectors.

        Args:
            film (Material): The film material.
            angle (float): The twist angle in degrees.
            merged_material (Material): The merged material.

        Returns:
            Material: The final material with the orthogonal supercell.
        """
        p, q = calculate_moire_periodicity(film.lattice.b, film.lattice.a, angle)
        new_atoms = filter_by_box(
            merged_material,
            [0, 0, 0],
            [p, q, merged_material.lattice.c],
            use_cartesian_coordinates=True,
        )

        new_lattice = new_atoms.lattice
        new_lattice.a = q
        new_lattice.b = p

        final_material = new_atoms.clone()
        final_material.set_new_lattice_vectors([p, 0, 0], [0, q, 0], [0, 0, new_lattice.c])

        return final_material

    @staticmethod
    def _calculate_supercell_size(angle: float, mod: float):
        """
        Reduces any given twist angle to the equivalent small angle.

        Parameters:
        angle (float): Twist angle in degrees.
        mod (float): The modulo value of the symmetry.

        Returns:
        float: Reduced twist angle within
        """
        mod_angle = angle % mod

        if mod_angle > mod / 2:
            mod_angle = mod - mod_angle

        n = np.round(1 / np.sin(np.deg2rad(mod_angle)))
        m = 2 * n
        return n, m

    def _update_material_name(self, material: Material, configuration: TwistedInterfaceConfiguration) -> Material:
        film_formula = get_chemical_formula(configuration.film_configuration.bulk)
        substrate_formula = get_chemical_formula(configuration.substrate_configuration.bulk)
        film_miller = "".join(map(str, configuration.film_configuration.miller_indices))
        substrate_miller = "".join(map(str, configuration.substrate_configuration.miller_indices))
        twist_str = f"Twist {configuration.twist_angle:.2f} degrees"
        new_name = f"{film_formula}({film_miller})-{substrate_formula}({substrate_miller}), {twist_str}"
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


class NanoRibbonTwistedInterfaceConfiguration(BaseConfiguration):
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

    film: Material
    substrate: Material
    twist_angle: float = 0
    ribbon_width: int = 1
    ribbon_length: int = 1
    distance_z: float = 3.0
    vacuum_x: int = 2
    vacuum_y: int = 2

    @property
    def _json(self):
        return {
            "type": self.get_cls_name(),
            "film": self.film.to_json(),
            "substrate": self.substrate.to_json(),
            "twist_angle": self.twist_angle,
            "ribbon_width": self.ribbon_width,
            "ribbon_length": self.ribbon_length,
            "distance_z": self.distance_z,
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
        top_ribbon = translate_by_vector(top_ribbon, translation_vector)
        merged_material = merge_materials([bottom_ribbon, top_ribbon])

        return [merged_material]

    def _update_material_name(
        self, material: Material, configuration: NanoRibbonTwistedInterfaceConfiguration
    ) -> Material:
        material.name = f"Twisted Nanoribbon Interface ({configuration.twist_angle:.2f}°)"
        return material


class CommensurateSuperCellTwistedInterfaceConfiguration(BaseConfiguration):
    film: Material
    substrate: Material
    target_angle: float = Field(..., description="Target twist angle in degrees")
    max_supercell_size: int = Field(10, description="Maximum supercell size to consider")
    distance_z: float = Field(3.0, description="Vertical distance between layers in Angstroms")

    @property
    def _json(self):
        return {
            "type": self.get_cls_name(),
            "film": self.film.to_json(),
            "substrate": self.substrate.to_json(),
            "target_angle": self.target_angle,
            "max_supercell_size": self.max_supercell_size,
            "distance_z": self.distance_z,
        }


class CommensurateSuperCellTwistedInterfaceBuilder(BaseBuilder):
    _ConfigurationType = CommensurateSuperCellTwistedInterfaceConfiguration

    def _generate(self, configuration: CommensurateSuperCellTwistedInterfaceConfiguration) -> List[Material]:
        n, m, p, q, actual_angle = self._find_commensurate_supercell(
            configuration.film, configuration.target_angle, configuration.max_supercell_size
        )

        # Create bottom supercell
        bottom_supercell = create_supercell(configuration.substrate, [[n, -m, 0], [m, n, 0], [0, 0, 1]])

        # Create top supercell and rotate it
        top_supercell = create_supercell(configuration.film, [[p, -q, 0], [q, p, 0], [0, 0, 1]])
        top_supercell = rotate_material(top_supercell, [0, 0, 1], actual_angle, wrap=False)

        # Translate top supercell
        translation_vector = [0, 0, configuration.distance_z]
        top_supercell = translate_by_vector(top_supercell, translation_vector)

        # Merge supercells
        final_material = merge_materials([bottom_supercell, top_supercell])

        return [final_material]

    def _find_commensurate_supercell(
        self, material: Material, target_angle: float, max_size: int
    ) -> Tuple[int, int, int, int, float]:
        best_angle_diff = float("inf")
        best_params = (1, 0, 1, 0, 0)

        for n in range(1, max_size + 1):
            for m in range(max_size + 1):
                for p in range(1, max_size + 1):
                    for q in range(max_size + 1):
                        angle = np.degrees(np.arctan2(n * q - m * p, n * p + m * q))
                        if abs(angle - target_angle) < best_angle_diff:
                            best_angle_diff = abs(angle - target_angle)
                            best_params = (n, m, p, q, angle)

        return best_params

    def _update_material_name(
        self, material: Material, configuration: CommensurateSuperCellTwistedInterfaceConfiguration
    ) -> Material:
        material.name = f"Commensurate Twisted Interface : {material.metadata['build']['actual_angle']:.2f} degrees"
        return material
