from typing import Any, Optional, Union, Type

from mat3ra.made.lattice import Lattice
from mat3ra.made.material import Material
from . import NanoTapeBuilderParameters
from .configuration import NanoTapeConfiguration
from .....build_components import MaterialWithBuildMetadata, TypeConfiguration
from .....build_components.entities.reusable.one_dimensional.crystal_lattice_lines.edge_types import EdgeTypesEnum
from .....build_components.entities.reusable.one_dimensional.crystal_lattice_lines_unique_repeated.builder import (
    CrystalLatticeLinesRepeatedBuilder,
)
from .....build_components.entities.reusable.one_dimensional.crystal_lattice_lines_unique_repeated import (
    CrystalLatticeLinesUniqueRepeatedConfiguration,
)
from .....build_components.operations.core.combinations.stack.builder import StackNComponentsBuilder
from .....modify import wrap_to_unit_cell


class NanoTapeBuilder(StackNComponentsBuilder):
    _ConfigurationType: Type[NanoTapeConfiguration] = NanoTapeConfiguration
    _GeneratedItemType: Type[Material] = Material
    _BuildParametersType: Type[NanoTapeBuilderParameters] = NanoTapeBuilderParameters
    _DefaultBuildParameters: NanoTapeBuilderParameters = NanoTapeBuilderParameters(
        use_rectangular_lattice=True,
    )

    @property
    def stack_component_types_conversion_map(self):
        return {
            **super().stack_component_types_conversion_map,
            CrystalLatticeLinesUniqueRepeatedConfiguration: CrystalLatticeLinesRepeatedBuilder,
        }

    def _make_rectangular_lattice(self, item: Material) -> Material:
        lattice_vectors = item.lattice.vector_arrays
        new_lattice_vectors = [
            [lattice_vectors[0][0], 0.0, 0.0],
            [0.0, lattice_vectors[1][1], 0.0],
            [0.0, 0.0, lattice_vectors[2][2]],
        ]
        new_lattice = Lattice.from_vectors_array(vectors=new_lattice_vectors)
        item.set_lattice(new_lattice)
        return item

    def _get_edge_type_from_miller_indices(self, miller_indices_2d: tuple) -> str:
        if miller_indices_2d == (1, 1):
            return EdgeTypesEnum.armchair.value.capitalize()
        elif miller_indices_2d == (0, 1):
            return EdgeTypesEnum.zigzag.value.capitalize()
        else:
            miller_str = f"{miller_indices_2d[0]}{miller_indices_2d[1]}"
            return f"({miller_str})"

    def _update_material_name_with_edge_type(
        self,
        material: Union[Material, MaterialWithBuildMetadata],
        crystal_name: str,
        miller_indices_2d: tuple,
        structure_type: str,
    ) -> Material:
        edge_type = self._get_edge_type_from_miller_indices(miller_indices_2d)
        miller_str = f"{miller_indices_2d[0]}{miller_indices_2d[1]}"
        material.name = f"{crystal_name} - {edge_type} {structure_type} ({miller_str})"
        return material

    def _update_material_name(
        self, material: Union[Material, MaterialWithBuildMetadata], configuration: Any
    ) -> Material:
        if isinstance(configuration, NanoTapeConfiguration):
            lattice_lines = configuration.lattice_lines
            material = self._update_material_name_with_edge_type(
                material, lattice_lines.crystal.name, lattice_lines.miller_indices_2d, "Nanotape"
            )
            return material
        return super()._update_material_name(material, configuration)

    def _post_process(
        self,
        item: Material,
        post_process_parameters: Optional[Any] = None,
        configuration: Optional[TypeConfiguration] = None,
    ) -> Material:
        item = super()._post_process(item, post_process_parameters, configuration)
        params = self.build_parameters or self._DefaultBuildParameters
        if params.use_rectangular_lattice:
            item = self._make_rectangular_lattice(item)
        return wrap_to_unit_cell(item)
