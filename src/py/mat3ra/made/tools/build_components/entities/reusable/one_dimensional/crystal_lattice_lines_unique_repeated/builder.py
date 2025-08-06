from typing import Type

from mat3ra.made.material import Material

from ......analyze.lattice_lines import CrystalLatticeLinesMaterialAnalyzer
from ......modify import wrap_to_unit_cell
from ......operations.core.unary import supercell, translate
from ..crystal_lattice_lines.builder import CrystalLatticeLinesBuilder
from .configuration import CrystalLatticeLinesUniqueRepeatedConfiguration


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
