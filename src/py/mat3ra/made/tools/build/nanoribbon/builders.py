from typing import List, Optional, Any, Type
import numpy as np
from mat3ra.made.tools.build.slab.utils import get_orthogonal_c_slab

from mat3ra.made.material import Material
from mat3ra.made.tools.build import BaseSingleBuilder
from ..stack.builders import Stack2ComponentsBuilder
from ...analyze.lattice_lines import CrystalLatticeLinesAnalyzer
from ...operations.core.unary import supercell, translate
from ...modify import wrap_to_unit_cell, translate_to_z_level, filter_by_rectangle_projection

from .configuration import (
    CrystalLatticeLinesConfiguration,
    CrystalLatticeLinesUniqueRepeatedConfiguration,
    NanoribbonConfiguration,
)


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

    _ConfigurationType: Type[CrystalLatticeLinesUniqueRepeatedConfiguration] = (
        CrystalLatticeLinesUniqueRepeatedConfiguration
    )

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


class NanoribbonBuilder(Stack2ComponentsBuilder):
    """
    Builder for creating nanoribbons from monolayer materials.
    Uses the new configuration structure with crystal lattice lines and vacuum stacking.
    """

    _ConfigurationType = NanoribbonConfiguration
    _GeneratedItemType = Material

    def _configuration_to_material(self, configuration_or_material: Any) -> Material:
        """Convert configuration objects to materials using dedicated builders."""
        if isinstance(configuration_or_material, CrystalLatticeLinesUniqueRepeatedConfiguration):
            builder = CrystalLatticeLinesRepeatedBuilder()
            return builder.get_material(configuration_or_material)
        # Fall back to parent method for other types (VacuumConfiguration, etc.)
        return super()._configuration_to_material(configuration_or_material)

    def _generate(self, configuration: NanoribbonConfiguration) -> Material:
        # Generate the basic stacked material (lattice lines + vacuum)
        stacked_material = super()._generate(configuration)

        # Apply nanoribbon-specific transformations
        nanoribbon = self._create_nanoribbon_from_stacked_material(stacked_material, configuration)

        return nanoribbon

    def _create_nanoribbon_from_stacked_material(
        self, stacked_material: Material, configuration: NanoribbonConfiguration
    ) -> Material:
        """
        Create the final nanoribbon by applying width/length repetitions and vacuum.
        """
        material = stacked_material.clone()

        # Calculate dimensions for the nanoribbon
        width_repetitions = configuration.width
        length_repetitions = configuration.length
        vacuum_width = configuration.vacuum_width
        vacuum_length = configuration.vacuum_length

        # Create repetitions in width and length directions
        # This creates the nanoribbon with the specified dimensions
        total_width_repetitions = width_repetitions + vacuum_width
        total_length_repetitions = length_repetitions + vacuum_length

        # Create a large supercell
        supercell_matrix = [[total_length_repetitions, 0, 0], [0, total_width_repetitions, 0], [0, 0, 1]]
        large_supercell = supercell(material, supercell_matrix)

        # Calculate the actual nanoribbon dimensions in cartesian coordinates
        lattice_vectors = large_supercell.lattice.vector_arrays
        length_cartesian = length_repetitions * lattice_vectors[0][0] / total_length_repetitions
        width_cartesian = width_repetitions * lattice_vectors[1][1] / total_width_repetitions
        height_cartesian = lattice_vectors[2][2]

        # TODO: do this via stacking vacuum on x and then on y
        # Update lattice to include vacuum regions
        new_lattice_vectors = [
            [length_cartesian + (vacuum_length * lattice_vectors[0][0] / total_length_repetitions), 0.0, 0.0],
            [0.0, width_cartesian + (vacuum_width * lattice_vectors[1][1] / total_width_repetitions), 0.0],
            [0.0, 0.0, height_cartesian],
        ]

        from mat3ra.made.lattice import Lattice

        new_lattice = Lattice.from_vectors_array(vectors=new_lattice_vectors)
        large_supercell.set_lattice(new_lattice)

        return large_supercell

    def _update_material_name(self, material: Material, configuration: NanoribbonConfiguration) -> Material:
        """Update material name to reflect nanoribbon properties."""
        lattice_lines = configuration.lattice_lines
        miller_str = f"{configuration.miller_indices_uv[0]}{configuration.miller_indices_uv[1]}"

        # Convert (u,v) to edge type for naming
        if configuration.miller_indices_uv == (1, 1):
            edge_type = "Armchair"
        elif configuration.miller_indices_uv == (0, 1):
            edge_type = "Zigzag"
        else:
            edge_type = f"({miller_str})"

        material.name = f"{lattice_lines.crystal.name} - {edge_type} Nanoribbon ({miller_str})"
        return material

    def _post_process(self, item: Material, post_process_parameters: Optional[Any] = None) -> Material:
        """Post-process the nanoribbon material."""
        item = super()._post_process(item, post_process_parameters)
        return wrap_to_unit_cell(item)
