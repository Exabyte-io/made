from typing import Any, List, Optional, Type

from mat3ra.code.entity import InMemoryEntityPydantic
from pydantic import BaseModel

from mat3ra.made.material import Material
from .configuration import (
    InterfaceConfiguration,
    NanoRibbonTwistedInterfaceConfiguration,
)
from ..nanoribbon import create_nanoribbon
from ..slab.builders import SlabStrainedSupercellBuilder
from ..slab.builders import SlabWithGapBuilder
from ..slab.configurations import SlabStrainedSupercellConfiguration
from ..slab.configurations import SlabStrainedSupercellWithGapConfiguration
from ..stack.builders import StackNComponentsBuilder
from ..stack.configuration import StackConfiguration
from ..utils import merge_materials
from ...analyze.other import get_chemical_formula
from ...build import BaseBuilder
from ...convert.utils import InterfacePartsEnum
from ...modify import (
    rotate,
    translate_by_vector,
    add_vacuum_sides,
    wrap_to_unit_cell,
)
from ...utils import AXIS_TO_INDEX_MAP


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
        if isinstance(configuration_or_material, SlabStrainedSupercellWithGapConfiguration):
            builder = SlabWithGapBuilder()
            return builder.get_material(configuration_or_material)
        elif isinstance(configuration_or_material, SlabStrainedSupercellConfiguration):
            builder = SlabStrainedSupercellBuilder()
            return builder.get_material(configuration_or_material)
        return super()._configuration_to_material(configuration_or_material)

    def _generate(self, configuration: InterfaceConfiguration) -> Material:
        film_material = self._configuration_to_material(configuration.film_configuration)
        substrate_material = self._configuration_to_material(configuration.substrate_configuration)

        film_material.set_labels_from_value(InterfacePartsEnum.FILM)
        substrate_material.set_labels_from_value(InterfacePartsEnum.SUBSTRATE)

        stacking_axis = AXIS_TO_INDEX_MAP[configuration.direction.value]
        other_axes = [i for i in range(3) if i != stacking_axis]
        translation_vector = [0.0, 0.0, 0.0]
        translation_vector[other_axes[0]] = configuration.xy_shift[0]
        translation_vector[other_axes[1]] = configuration.xy_shift[1]

        film_material = translate_by_vector(film_material, translation_vector, use_cartesian_coordinates=True)
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
#                       Twisted Interface Builders                                     #
########################################################################################


class NanoRibbonTwistedInterfaceBuilder(BaseBuilder):
    _GeneratedItemType = Material
    _ConfigurationType: type(  # type: ignore
        NanoRibbonTwistedInterfaceConfiguration
    ) = NanoRibbonTwistedInterfaceConfiguration  # type: ignore

    def _generate(self, configuration: _ConfigurationType) -> List[Material]:
        bottom_ribbon = create_nanoribbon(
            material=configuration.substrate,
            width=configuration.ribbon_width,
            length=configuration.ribbon_length,
        )
        top_ribbon = create_nanoribbon(
            material=configuration.film,
            width=configuration.ribbon_width,
            length=configuration.ribbon_length,
        )
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


class CommensurateLatticeInterfaceBuilderParameters(BaseModel):
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
