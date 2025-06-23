from typing import Any, List, Optional, Type

import numpy as np
from mat3ra.code.entity import InMemoryEntityPydantic
from pydantic import BaseModel
from pymatgen.analysis.interfaces.coherent_interfaces import (
    CoherentInterfaceBuilder,
    ZSLGenerator,
)

from mat3ra.made.material import Material

from .configuration import (
    InterfaceConfiguration,
    NanoRibbonTwistedInterfaceConfiguration,
)
from .enums import StrainModes
from .termination_pair import TerminationPair, safely_select_termination_pair
from .utils import interface_patch_with_mean_abs_strain, remove_duplicate_interfaces
from ..mixins import (
    ConvertGeneratedItemsPymatgenStructureMixin,
)
from ..nanoribbon import NanoribbonConfiguration, create_nanoribbon
from ..slab.builders import SlabStrainedSupercellBuilder
from ..slab.configuration import SlabStrainedSupercellConfiguration
from ..stack.builders import StackNComponentsBuilder
from ..stack.configuration import StackConfiguration
from ..utils import merge_materials
from ...analyze.other import get_chemical_formula
from ...build import BaseBuilder, BaseBuilderParameters
from ...convert import to_pymatgen, PymatgenInterface
from ...modify import (
    translate_to_z_level,
    rotate,
    translate_by_vector,
    add_vacuum_sides,
    wrap_to_unit_cell,
)


class InterfaceBuilderParameters(InMemoryEntityPydantic):
    pass


class InterfaceBuilder(StackNComponentsBuilder):
    """
    Creates matching interface between substrate and film by straining the film to match the substrate.
    """

    _ConfigurationType = InterfaceConfiguration
    _BuildParametersType = InterfaceBuilderParameters
    _GeneratedItemType: Type[Material] = Material

    def _configuration_to_material(self, configuration_or_material: Any) -> Material:
        if isinstance(configuration_or_material, SlabStrainedSupercellConfiguration):
            builder = SlabStrainedSupercellBuilder()
            return builder.get_material(configuration_or_material)
        return super()._configuration_to_material(configuration_or_material)

    def _generate(self, configuration: InterfaceConfiguration) -> Material:
        film_material = self._configuration_to_material(configuration.film_configuration)
        substrate_material = self._configuration_to_material(configuration.substrate_configuration)

        film_material.set_labels_from_value(0)
        substrate_material.set_labels_from_value(1)

        film_material = translate_by_vector(film_material, configuration.xy_shift + [0], use_cartesian_coordinates=True)
        stack_configuration = StackConfiguration(
            stack_components=[substrate_material, film_material, configuration.vacuum_configuration],
            direction=configuration.direction,
        )

        interface = super()._generate(stack_configuration)

        wrapped_interface = wrap_to_unit_cell(interface)
        return wrapped_interface

    def _update_material_name(self, material: Material, configuration: InterfaceConfiguration) -> Material:
        film_formula = get_chemical_formula(configuration.film_configuration.atomic_layers.crystal)
        substrate_formula = get_chemical_formula(configuration.substrate_configuration.atomic_layers.crystal)
        film_miller_indices = "".join([str(i) for i in configuration.film_configuration.atomic_layers.miller_indices])
        substrate_miller_indices = "".join(
            [str(i) for i in configuration.substrate_configuration.atomic_layers.miller_indices]
        )
        name = f"{film_formula}({film_miller_indices})-{substrate_formula}({substrate_miller_indices}), Interface"
        material.name = name
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
