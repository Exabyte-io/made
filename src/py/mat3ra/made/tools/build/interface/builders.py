from typing import Type, Union

import numpy as np
from mat3ra.code.entity import InMemoryEntityPydantic

from mat3ra.made.material import Material
from .configuration import InterfaceConfiguration
from .. import MaterialWithBuildMetadata
from ..slab.builders import SlabStrainedSupercellBuilder
from ..slab.configurations import SlabStrainedSupercellConfiguration
from ..stack.builders import StackNComponentsBuilder
from ..stack.configuration import StackConfiguration
from ...analyze import BaseMaterialAnalyzer
from ...analyze.lattice import get_material_with_primitive_lattice
from ...convert.utils import InterfacePartsEnum
from ...modify import (
    translate_by_vector,
    translate_to_z_level,
    wrap_to_unit_cell,
)
from ...operations.core.unary import supercell
from ...operations.core.utils import should_skip_stacking
from ....utils import AXIS_TO_INDEX_MAP


class InterfaceBuilderParameters(InMemoryEntityPydantic):
    make_primitive: bool = True


class InterfaceBuilder(StackNComponentsBuilder):
    """
    Creates an interface by straining the film to match the substrate.
    """

    _ConfigurationType = InterfaceConfiguration
    _BuildParametersType = InterfaceBuilderParameters
    _DefaultBuildParameters = InterfaceBuilderParameters()
    _GeneratedItemType: Type[Material] = Material

    @property
    def stack_component_types_conversion_map(self):
        return {
            **super().stack_component_types_conversion_map,
            SlabStrainedSupercellConfiguration: SlabStrainedSupercellBuilder,
        }

    def _generate(self, configuration: InterfaceConfiguration) -> MaterialWithBuildMetadata:
        film_material = self._stack_component_to_material(configuration.film_configuration, configuration)
        substrate_material = self._stack_component_to_material(configuration.substrate_configuration, configuration)

        film_material.set_labels_from_value(InterfacePartsEnum.FILM.value)
        substrate_material.set_labels_from_value(InterfacePartsEnum.SUBSTRATE.value)

        # Apply xy shift to the film material
        stacking_axis = AXIS_TO_INDEX_MAP[configuration.direction.value]
        other_axes = [i for i in range(3) if i != stacking_axis]
        translation_vector = [0.0, 0.0, 0.0]
        translation_vector[other_axes[0]] = configuration.xy_shift[0]
        translation_vector[other_axes[1]] = configuration.xy_shift[1]
        film_material = translate_by_vector(film_material, translation_vector, use_cartesian_coordinates=True)

        try:
            should_skip_stacking(substrate_material, film_material, stacking_axis)
        except ValueError:
            film_material = supercell(film_material, [[0, 1, 0], [1, 0, 0], [0, 0, 1]])
            print(
                "Switched in-plane lattice vectors of the film material to match substrate"
                + "Direct stacking was not possible."
            )
            should_skip_stacking(substrate_material, film_material, stacking_axis)

        stack_configuration = StackConfiguration(
            stack_components=[substrate_material, film_material, configuration.vacuum_configuration],
            gaps=configuration.gaps,
            direction=configuration.direction,
        )
        interface = super()._generate(stack_configuration)
        wrapped_interface = wrap_to_unit_cell(interface)
        return wrapped_interface

    def _post_process(
        self, material: MaterialWithBuildMetadata, configuration: InterfaceConfiguration
    ) -> MaterialWithBuildMetadata:
        if self.build_parameters.make_primitive:
            # TODO: check that this doesn't warp material or flip it -- otherwise raise and skip
            primitive_material = get_material_with_primitive_lattice(material, return_original_if_not_reduced=True)

            if primitive_material != material:
                primitive_material = self._preserve_interface_labels_after_primitive_conversion(
                    original_material=material, primitive_material=primitive_material
                )

            material = primitive_material

        return super()._post_process(material, configuration)

    def get_base_name_from_configuration(self, first_configuration, second_configuration) -> str:
        first_component_formula = BaseMaterialAnalyzer(material=first_configuration.atomic_layers.crystal).formula
        second_component_formula = BaseMaterialAnalyzer(material=second_configuration.atomic_layers.crystal).formula
        first_component_miller = first_configuration.atomic_layers.miller_indices_as_string
        second_component_miller = second_configuration.atomic_layers.miller_indices_as_string
        return f"{first_component_formula}{first_component_miller}-{second_component_formula}{second_component_miller}"

    def get_name_suffix(self, configuration: InterfaceConfiguration) -> str:
        strain = configuration.von_mises_strain_percentage
        if strain == 0:
            return "Interface"
        return f"Interface, Strain {strain:.3f}pct"

    def _update_material_name(
        self, material: Union[Material, MaterialWithBuildMetadata], configuration: InterfaceConfiguration
    ) -> MaterialWithBuildMetadata:
        base_name = self.get_base_name_from_configuration(
            configuration.film_configuration, configuration.substrate_configuration
        )
        name_suffix = self.get_name_suffix(configuration) or self.name_suffix
        name = f"{base_name}, {name_suffix}"
        material.name = name
        return material

    def _preserve_interface_labels_after_primitive_conversion(
        self, original_material: MaterialWithBuildMetadata, primitive_material: MaterialWithBuildMetadata
    ) -> MaterialWithBuildMetadata:
        if not original_material.basis.labels.values:
            return primitive_material

        labeled_primitive = primitive_material.clone()
        interface_z = self._find_interface_z_level(original_material)
        primitive_labels = self._assign_labels_by_z_level(labeled_primitive, interface_z)
        labeled_primitive.set_labels_from_list(primitive_labels)
        return labeled_primitive

    def _find_interface_z_level(self, material: MaterialWithBuildMetadata) -> float:
        normalized_material = translate_to_z_level(material, z_level="bottom")
        normalized_material.to_cartesian()
        z_coords = np.array(normalized_material.coordinates_array)[:, 2]
        labels = material.basis.labels.values

        substrate_z = z_coords[np.array(labels) == int(InterfacePartsEnum.SUBSTRATE.value)]
        film_z = z_coords[np.array(labels) == int(InterfacePartsEnum.FILM.value)]

        if len(substrate_z) == 0 or len(film_z) == 0:
            return (np.min(z_coords) + np.max(z_coords)) / 2

        max_substrate_z = np.max(substrate_z)
        min_film_z = np.min(film_z)
        return (max_substrate_z + min_film_z) / 2

    def _assign_labels_by_z_level(self, material: MaterialWithBuildMetadata, interface_z: float) -> list:
        normalized_material = translate_to_z_level(material, z_level="bottom")
        normalized_material.to_cartesian()
        z_coords = np.array(normalized_material.coordinates_array)[:, 2]
        substrate_label = int(InterfacePartsEnum.SUBSTRATE.value)
        film_label = int(InterfacePartsEnum.FILM.value)
        return [substrate_label if z < interface_z else film_label for z in z_coords]
