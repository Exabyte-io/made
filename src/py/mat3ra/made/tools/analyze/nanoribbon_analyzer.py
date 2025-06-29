from typing import Optional

from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from mat3ra.made.tools.analyze.base_stack_analyzer import BaseStackAnalyzer
from mat3ra.made.tools.build.nanoribbon.nano_tape_configuration import NanoTapeConfiguration
from mat3ra.made.tools.build.nanoribbon.nanoribbon_configuration import NanoribbonConfiguration
from mat3ra.made.tools.build.slab.entities import Termination


class NanoribbonAnalyzer(BaseStackAnalyzer[NanoTapeConfiguration]):
    """
    Analyzer for creating nanoribbon configurations.
    Creates NanoribbonConfiguration with nanotape and vacuum.
    """

    @property
    def nanotape_analyzer(self):
        from mat3ra.made.tools.analyze.nanotape_analyzer import NanoTapeAnalyzer

        return NanoTapeAnalyzer(material=self.material, miller_indices_uv=self.miller_indices_uv)

    def get_component_configuration(
        self,
        width: int,
        length: int,
        termination: Optional[Termination] = None,
        vacuum_width: float = 10.0,
    ) -> NanoTapeConfiguration:
        """Get the nanotape configuration."""
        return self.nanotape_analyzer.get_configuration(
            width=width,
            length=length,
            termination=termination,
            vacuum_width=vacuum_width,
        )

    def get_component_builder(self):
        from mat3ra.made.tools.build.nanoribbon.builders import NanoTapeBuilder

        return NanoTapeBuilder()

    def get_stack_direction(self) -> AxisEnum:
        return AxisEnum.x

    def get_vacuum_direction(self) -> AxisEnum:
        return AxisEnum.x

    def get_configuration(
        self,
        width: int,
        length: int,
        termination: Optional[Termination] = None,
        vacuum_width: float = 10.0,
        vacuum_length: float = 10.0,
    ) -> NanoribbonConfiguration:
        """
        Create a nanoribbon configuration with nanotape and vacuum.
        Args:
            width: The width of the nanoribbon in number of unit cells.
            length: The length of the nanoribbon in number of unit cells.
            termination: The termination to use for the nanoribbon.
            vacuum_width: The width of the vacuum region in Angstroms (cartesian).
            vacuum_length: The length of the vacuum region in Angstroms (cartesian).
        """
        stack_config = super().get_configuration(
            vacuum_size=vacuum_length,
            width=width,
            length=length,
            termination=termination,
            vacuum_width=vacuum_width,
        )
        return NanoribbonConfiguration(**stack_config.model_dump())
