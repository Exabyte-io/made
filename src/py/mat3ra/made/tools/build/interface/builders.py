from typing import Any, Type

from mat3ra.code.entity import InMemoryEntityPydantic

from mat3ra.made.material import Material
from .configuration import (
    InterfaceConfiguration,
)
from ..slab.builders import SlabStrainedSupercellBuilder
from ..slab.builders import SlabWithGapBuilder
from ..slab.configurations import SlabStrainedSupercellConfiguration
from ..slab.configurations import SlabStrainedSupercellWithGapConfiguration
from ..stack.builders import StackNComponentsBuilder
from ..stack.configuration import StackConfiguration
from ...analyze.other import get_chemical_formula
from ...convert.utils import InterfacePartsEnum
from ...modify import (
    translate_by_vector,
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
