from typing import Optional, Tuple
from mat3ra.made.material import Material
from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from mat3ra.made.tools.analyze.lattice import LatticeMaterialAnalyzer
from mat3ra.made.tools.analyze.lattice_lines import CrystalLatticeLinesAnalyzer
from mat3ra.made.tools.build.vacuum.configuration import VacuumConfiguration
from mat3ra.made.tools.build.slab.entities import Termination
from mat3ra.made.tools.build.nanoribbon.configuration import CrystalLatticeLinesUniqueRepeatedConfiguration, NanoTapeConfiguration

class NanoTapeAnalyzer(LatticeMaterialAnalyzer):
    """
    Analyzer for creating nanotape configurations.
    Creates NanoTapeConfiguration with repeated crystal lattice lines and vacuum.
    """
    miller_indices_uv: Tuple[int, int]

    def __init__(self, monolayer: Material, miller_indices_uv: Tuple[int, int]):
        super().__init__(material=monolayer, miller_indices_uv=miller_indices_uv)

    @property
    def lattice_lines_analyzer(self):
        return CrystalLatticeLinesAnalyzer(self.material, self.miller_indices_uv)

    def get_configuration(
        self,
        width: int,
        length: int,
        termination: Optional[Termination] = None,
        vacuum_width: float = 10.0,
    ) -> NanoTapeConfiguration:
        """
        Create a nanotape configuration with repeated lattice lines and vacuum.
        Args:
            width: The width of the nanotape in number of unit cells.
            length: The length of the nanotape in number of unit cells.
            termination: The termination to use for the nanotape.
            vacuum_width: The width of the vacuum region in Angstroms (cartesian).
        """
        if termination is None:
            termination = self.lattice_lines_analyzer.default_termination
        lattice_lines_config = CrystalLatticeLinesUniqueRepeatedConfiguration(
            crystal=self.material,
            miller_indices_uv=self.miller_indices_uv,
            termination_top=termination,
            number_of_repetitions_width=width,
            number_of_repetitions_length=length,
        )
        from mat3ra.made.tools.build.nanoribbon.builders import CrystalLatticeLinesRepeatedBuilder
        lattice_lines_builder = CrystalLatticeLinesRepeatedBuilder()
        lattice_lines_material = lattice_lines_builder.get_material(lattice_lines_config)
        vacuum_config = VacuumConfiguration(
            size=vacuum_width,
            crystal=lattice_lines_material,
            direction=AxisEnum.y,
        )
        return NanoTapeConfiguration(
            stack_components=[lattice_lines_config, vacuum_config],
            direction=AxisEnum.y,
        ) 