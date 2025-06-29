from typing import Optional, Any, Type

from .enums import EdgeTypes
from pydantic import Field

from mat3ra.made.lattice import Lattice
from mat3ra.made.material import Material
from mat3ra.made.tools.build.slab.utils import get_orthogonal_c_slab
from .crystal_lattice_lines_configuration import (
    CrystalLatticeLinesConfiguration,
    CrystalLatticeLinesUniqueRepeatedConfiguration,
)
from .nano_tape_configuration import NanoTapeConfiguration
from .nanoribbon_configuration import NanoribbonConfiguration
from .. import BaseSingleBuilder, BaseBuilderParameters
from ..stack.builders import Stack2ComponentsBuilder
from ...analyze.lattice_lines import CrystalLatticeLinesAnalyzer
from ...modify import wrap_to_unit_cell, translate_to_z_level
from ...operations.core.unary import supercell, translate


class CrystalLatticeLinesBuilder(BaseSingleBuilder):
    """
    Builder for creating a single crystal lattice line with termination.
    This is similar to CrystalLatticePlanesBuilder but for 1D lines.
    """

    _PostProcessParametersType: Any = None
    use_enforce_convention: bool = True

    def _generate(self, configuration: CrystalLatticeLinesConfiguration) -> Material:
        crystal_lattice_lines_analyzer = CrystalLatticeLinesAnalyzer(
            material=configuration.crystal, miller_indices_uv=configuration.miller_indices_uv
        )
        miller_supercell_matrix = crystal_lattice_lines_analyzer.miller_supercell_matrix
        miller_supercell_material = supercell(configuration.crystal, miller_supercell_matrix)
        # Lattice returned with vector B being treated as the direction to the surface, needs to be swaped with C
        rotated_material = supercell(miller_supercell_material, [[1, 0, 0], [0, 0, 1], [0, 1, 0]])
        orthogonal_material = get_orthogonal_c_slab(rotated_material)
        return orthogonal_material

    def _enforce_convention(self, material: Material) -> Material:
        if not self.use_enforce_convention:
            return material
        return translate_to_z_level(material, "bottom")

    def _post_process(self, item: Material, post_process_parameters: Optional[_PostProcessParametersType]) -> Material:
        item = super()._post_process(item, post_process_parameters)
        return self._enforce_convention(item)


class CrystalLatticeLinesRepeatedBuilder(CrystalLatticeLinesBuilder):
    """
    Builder for creating repeated crystal lattice lines with termination.
    This is similar to AtomicLayersUniqueRepeatedBuilder but for 1D lines.
    """

    _ConfigurationType: Type[
        CrystalLatticeLinesUniqueRepeatedConfiguration
    ] = CrystalLatticeLinesUniqueRepeatedConfiguration

    def _generate(self, configuration: CrystalLatticeLinesUniqueRepeatedConfiguration) -> Material:
        crystal_lattice_lines_material = super()._generate(configuration)

        crystal_lattice_lines_analyzer = CrystalLatticeLinesAnalyzer(
            material=configuration.crystal, miller_indices_uv=configuration.miller_indices_uv
        )
        translation_vector = crystal_lattice_lines_analyzer.get_translation_vector_for_termination_without_vacuum(
            configuration.termination_top
        )
        material_translated = translate(crystal_lattice_lines_material, translation_vector)
        material_translated_wrapped = wrap_to_unit_cell(material_translated)

        material_translated_wrapped_repeated = supercell(
            material_translated_wrapped,
            [
                [configuration.number_of_repetitions_length, 0, 0],
                [0, configuration.number_of_repetitions_width, 0],
                [0, 0, 1],
            ],
        )
        return material_translated_wrapped_repeated


class RectangularLatticeBuilderParameters(BaseBuilderParameters):
    use_rectangular_lattice: bool = Field(
        True, description="If True, set the XY lattice to be rectangular after stacking."
    )


class BaseRectangularLatticeBuilder(Stack2ComponentsBuilder):
    _BuilderParametersType = RectangularLatticeBuilderParameters
    _DefaultBuildParameters = RectangularLatticeBuilderParameters(use_rectangular_lattice=True)

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

    def _post_process(self, item: Material, post_process_parameters: Optional[Any] = None) -> Material:
        item = super()._post_process(item, post_process_parameters)
        # Use default parameters if build_parameters is None
        params = self.build_parameters or self._DefaultBuildParameters
        if params.use_rectangular_lattice:
            item = self._make_rectangular_lattice(item)
        return wrap_to_unit_cell(item)

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


class NanoTapeBuilderParameters(RectangularLatticeBuilderParameters):
    pass


class NanoribbonBuilderParameters(RectangularLatticeBuilderParameters):
    pass


class NanoTapeBuilder(BaseRectangularLatticeBuilder):
    _ConfigurationType = NanoTapeConfiguration
    _GeneratedItemType = Material
    _BuilderParametersType = NanoTapeBuilderParameters

    def _configuration_to_material(self, configuration_or_material: Any) -> Material:
        if isinstance(configuration_or_material, CrystalLatticeLinesUniqueRepeatedConfiguration):
            builder = CrystalLatticeLinesRepeatedBuilder()
            return builder.get_material(configuration_or_material)
        return super()._configuration_to_material(configuration_or_material)

    def _update_material_name(self, material: Material, configuration: NanoTapeConfiguration) -> Material:
        lattice_lines = configuration.lattice_lines
        material = self._update_material_name_with_edge_type(
            material, lattice_lines.crystal.name, lattice_lines.miller_indices_uv, "Nanotape"
        )
        return material


class NanoribbonBuilder(BaseRectangularLatticeBuilder):
    _ConfigurationType = NanoribbonConfiguration
    _GeneratedItemType = Material
    _BuilderParametersType = NanoribbonBuilderParameters

    def _configuration_to_material(self, configuration_or_material: Any) -> Material:
        if isinstance(configuration_or_material, NanoTapeConfiguration):
            builder = NanoTapeBuilder()
            return builder.get_material(configuration_or_material)
        return super()._configuration_to_material(configuration_or_material)

    def _update_material_name(self, material: Material, configuration: NanoribbonConfiguration) -> Material:
        nanotape = configuration.nanotape
        material = self._update_material_name_with_edge_type(
            material, nanotape.lattice_lines.crystal.name, nanotape.lattice_lines.miller_indices_uv, "Nanoribbon"
        )
        return material
