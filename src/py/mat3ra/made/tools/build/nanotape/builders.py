from typing import Any, Optional

from mat3ra.made.lattice import Lattice
from mat3ra.made.material import Material
from mat3ra.made.tools.build import BaseBuilderParameters
from pydantic import Field
from ..stack.builders import Stack2ComponentsBuilder
from ..lattice_lines.configuration import CrystalLatticeLinesUniqueRepeatedConfiguration, EdgeTypes
from ..lattice_lines.builders import CrystalLatticeLinesRepeatedBuilder
from .configuration import NanoTapeConfiguration
from ...modify import wrap_to_unit_cell


class NanoTapeBuilderParameters(BaseBuilderParameters):
    use_rectangular_lattice: bool = Field(
        True, description="If True, set the XY lattice to be rectangular after stacking."
    )


class NanoTapeBuilder(Stack2ComponentsBuilder):
    _ConfigurationType = "NanoTapeConfiguration"  # String type annotation to avoid circular import
    _GeneratedItemType = Material
    _BuilderParametersType = NanoTapeBuilderParameters
    _DefaultBuildParameters = NanoTapeBuilderParameters(
        use_rectangular_lattice=True,
    )

    def _configuration_to_material(self, configuration_or_material: Any) -> Material:
        if isinstance(configuration_or_material, NanoTapeConfiguration):
            return self.get_material(configuration_or_material)
        if isinstance(configuration_or_material, CrystalLatticeLinesUniqueRepeatedConfiguration):
            builder = CrystalLatticeLinesRepeatedBuilder()
            return builder.get_material(configuration_or_material)
        return super()._configuration_to_material(configuration_or_material)

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

    def _get_edge_type_from_miller_indices(self, miller_indices_uv: tuple) -> str:
        if miller_indices_uv == (1, 1):
            return EdgeTypes.armchair.value.capitalize()
        elif miller_indices_uv == (0, 1):
            return EdgeTypes.zigzag.value.capitalize()
        else:
            miller_str = f"{miller_indices_uv[0]}{miller_indices_uv[1]}"
            return f"({miller_str})"

    def _update_material_name_with_edge_type(
        self, material: Material, crystal_name: str, miller_indices_uv: tuple, structure_type: str
    ) -> Material:
        edge_type = self._get_edge_type_from_miller_indices(miller_indices_uv)
        miller_str = f"{miller_indices_uv[0]}{miller_indices_uv[1]}"
        material.name = f"{crystal_name} - {edge_type} {structure_type} ({miller_str})"
        return material

    def _update_material_name(self, material: Material, configuration: Any) -> Material:
        if isinstance(configuration, NanoTapeConfiguration):
            lattice_lines = configuration.lattice_lines
            material = self._update_material_name_with_edge_type(
                material, lattice_lines.crystal.name, lattice_lines.miller_indices_uv, "Nanotape"
            )
            return material
        return super()._update_material_name(material, configuration)

    def _post_process(self, item: Material, post_process_parameters: Optional[Any] = None) -> Material:
        item = super()._post_process(item, post_process_parameters)
        params = self.build_parameters or self._DefaultBuildParameters
        if params.use_rectangular_lattice:
            item = self._make_rectangular_lattice(item)
        return wrap_to_unit_cell(item)
