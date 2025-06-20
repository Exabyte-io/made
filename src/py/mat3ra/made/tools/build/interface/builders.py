from typing import Any, List, Optional, Type

import numpy as np
from mat3ra.code.entity import InMemoryEntityPydantic
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
from ..metadata import MaterialMetadata
from ..mixins import (
    ConvertGeneratedItemsPymatgenStructureMixin,
)
from ..nanoribbon import NanoribbonConfiguration, create_nanoribbon
from ..slab.builders import SlabStrainedSupercellBuilder
from ..slab.configuration import SlabStrainedSupercellConfiguration
from ..stack.builders import StackBuilder2Components
from ..stack.configuration import StackConfiguration
from ..supercell import create_supercell
from ..utils import merge_materials
from ...build import BaseBuilder, BaseBuilderParameters
from ...convert import to_pymatgen, PymatgenInterface
from ...modify import (
    translate_to_z_level,
    rotate,
    translate_by_vector,
    add_vacuum_sides,
    add_vacuum,
    wrap_to_unit_cell,
)


class InterfaceBuilderParameters(InMemoryEntityPydantic):
    pass


class InterfaceBuilder(StackBuilder2Components):
    """
    Creates matching interface between substrate and film by straining the film to match the substrate.
    """

    _ConfigurationType = InterfaceConfiguration
    _BuildParametersType = InterfaceBuilderParameters
    _GeneratedItemType: Type[Material] = Material

    def configuration_to_material(self, configuration_or_material: Any) -> Material:
        if isinstance(configuration_or_material, Material):
            return configuration_or_material

        if isinstance(configuration_or_material, SlabStrainedSupercellConfiguration):
            builder = SlabStrainedSupercellBuilder()
            return builder.get_material(configuration_or_material)

        return super().configuration_to_material(configuration_or_material)

    def _generate(self, configuration: InterfaceConfiguration) -> List[Material]:
        substrate_material = self.configuration_to_material(configuration.substrate_configuration)
        film_material = self.configuration_to_material(configuration.film_configuration)
        if configuration.xy_shift:
            film_material = translate_by_vector(
                film_material, configuration.xy_shift + [0], use_cartesian_coordinates=True
            )

        substrate_labeled = substrate_material.clone()
        substrate_labeled.set_labels([0] * len(substrate_labeled.basis.elements.values))

        film_labeled = film_material.clone()
        film_labeled.set_labels([1] * len(film_labeled.basis.elements.values))

        substrate_film_stack_config = StackConfiguration(
            stack_components=[substrate_labeled, film_labeled],
            direction=configuration.direction,
        )
        substrate_film_materials = super()._generate(substrate_film_stack_config)
        interface = substrate_film_materials[0]

        vacuum_config = configuration.vacuum_configuration
        if vacuum_config.size > 0:
            interface_vacuum_stack_config = StackConfiguration(
                stack_components=[interface, vacuum_config], direction=configuration.direction
            )
            interface_vacuum_materials = super()._generate(interface_vacuum_stack_config)
            interface = interface_vacuum_materials[0]

        wrapped_interface = wrap_to_unit_cell(interface)
        return [wrapped_interface]

    def _update_material_name(self, material: Material, configuration: InterfaceConfiguration) -> Material:
        """Update the material name to reflect that it's an interface."""
        material.name = configuration.name
        return material


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
            film_thickness=configuration.film_configuration.number_of_layers,
            substrate_thickness=configuration.substrate_configuration.number_of_layers,
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
            interface.metadata = {"actual_twist_angle": item.angle}
            if item.configuration.vacuum != 0:
                interface = add_vacuum(interface, item.configuration.vacuum)
            interfaces.append(interface)
        return interfaces

    def _update_material_metadata(self, material, configuration) -> Material:
        updated_material = super()._update_material_metadata(material, configuration)
        if "actual_twist_angle" in material.metadata:
            new_metadata = MaterialMetadata(**material.metadata)
            new_metadata.build.configuration.update(actual_twist_angle=material.metadata["actual_twist_angle"])
            updated_material.metadata = new_metadata.to_dict()
        return updated_material

    def _update_material_name(
        self, material: Material, configuration: NanoRibbonTwistedInterfaceConfiguration
    ) -> Material:
        material.name = f"Twisted Bilayer Interface ({material.metadata['actual_twist_angle']:.2f} degrees)"
        return material
