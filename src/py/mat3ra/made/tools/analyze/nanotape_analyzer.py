from typing import Optional

from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from mat3ra.made.tools.analyze.base_stack_analyzer import BaseStackAnalyzer
from mat3ra.made.tools.analyze.lattice_lines import CrystalLatticeLinesAnalyzer
from mat3ra.made.tools.build.nanoribbon.configuration import (
    CrystalLatticeLinesUniqueRepeatedConfiguration,
    NanoTapeConfiguration,
)
from mat3ra.made.tools.build.slab.entities import Termination


class NanoTapeAnalyzer(BaseStackAnalyzer[CrystalLatticeLinesUniqueRepeatedConfiguration]):
    """
    Analyzer for creating nanotape configurations.
    Creates NanoTapeConfiguration with repeated crystal lattice lines and vacuum.
    """

    @property
    def lattice_lines_analyzer(self):
        return CrystalLatticeLinesAnalyzer(material=self.material, miller_indices_uv=self.miller_indices_uv)

    def get_component_configuration(
        self,
        width: int,
        length: int,
        termination: Optional[Termination] = None,
    ) -> CrystalLatticeLinesUniqueRepeatedConfiguration:
        if termination is None:
            termination = self.lattice_lines_analyzer.default_termination
        return CrystalLatticeLinesUniqueRepeatedConfiguration(
            crystal=self.material,
            miller_indices_uv=self.miller_indices_uv,
            termination_top=termination,
            number_of_repetitions_width=width,
            number_of_repetitions_length=length,
        )

    def get_component_builder(self):
        from mat3ra.made.tools.build.nanoribbon.builders import CrystalLatticeLinesRepeatedBuilder

        return CrystalLatticeLinesRepeatedBuilder()

    def get_stack_direction(self) -> AxisEnum:
        return AxisEnum.y

    def get_vacuum_direction(self) -> AxisEnum:
        return AxisEnum.y

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
        stack_config = super().get_configuration(
            vacuum_size=vacuum_width,
            width=width,
            length=length,
            termination=termination,
        )
        return NanoTapeConfiguration(**stack_config.model_dump())
