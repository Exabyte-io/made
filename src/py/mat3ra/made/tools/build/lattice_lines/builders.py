from typing import Optional, Any, Type, Union

from mat3ra.made.material import Material
from mat3ra.made.tools.build.slab.utils import get_orthogonal_c_slab
from .configuration import (
    CrystalLatticeLinesConfiguration,
    CrystalLatticeLinesUniqueRepeatedConfiguration,
)
from .. import BaseSingleBuilder, MaterialWithBuildMetadata
from ...analyze.lattice_lines import CrystalLatticeLinesMaterialAnalyzer
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
        crystal_lattice_lines_analyzer = CrystalLatticeLinesMaterialAnalyzer(
            material=configuration.crystal, miller_indices_2d=configuration.miller_indices_2d
        )
        miller_supercell_matrix = crystal_lattice_lines_analyzer.miller_supercell_matrix
        miller_supercell_material = supercell(configuration.crystal, miller_supercell_matrix)
        # Lattice returned with vector B being treated as the direction to the surface, needs to be swaped with C
        rotated_material = supercell(miller_supercell_material, [[1, 0, 0], [0, 0, 1], [0, 1, 0]])
        orthogonal_material = get_orthogonal_c_slab(rotated_material)
        return orthogonal_material

    def _enforce_convention(self, material: Union[Material, MaterialWithBuildMetadata]) -> Material:
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

        crystal_lattice_lines_analyzer = CrystalLatticeLinesMaterialAnalyzer(
            material=configuration.crystal, miller_indices_2d=configuration.miller_indices_2d
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
