from typing import Optional, Tuple
from mat3ra.made.material import Material
from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from mat3ra.made.tools.analyze.lattice import LatticeMaterialAnalyzer
from mat3ra.made.tools.build.vacuum.configuration import VacuumConfiguration
from mat3ra.made.tools.build.slab.entities import Termination
from mat3ra.made.tools.build.nanoribbon.configuration import NanoTapeConfiguration, NanoribbonConfiguration
from mat3ra.made.tools.analyze.nanotape_analyzer import NanoTapeAnalyzer

class NanoribbonAnalyzer(LatticeMaterialAnalyzer):
    """
    Analyzer for creating nanoribbon configurations.
    Creates NanoribbonConfiguration with nanotape and vacuum.
    """
    miller_indices_uv: Tuple[int, int]

    def __init__(self, monolayer: Material, miller_indices_uv: Tuple[int, int]):
        super().__init__(material=monolayer, miller_indices_uv=miller_indices_uv)

    @property
    def nanotape_analyzer(self):
        return NanoTapeAnalyzer(self.material, self.miller_indices_uv)

    def get_configuration(
        self,
        width: int,
        length: int,
        termination: Optional[Termination] = None,
        vacuum_width: float = 10.0,
        vacuum_length: float = 0.0,
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
        nanotape_config = self.nanotape_analyzer.get_configuration(
            width=width,
            length=length,
            termination=termination,
            vacuum_width=vacuum_width,
        )
        
        # Always create vacuum configuration (VacuumBuilder handles zero case)
        from mat3ra.made.tools.build.nanoribbon.builders import NanoTapeBuilder
        nanotape_builder = NanoTapeBuilder()
        nanotape_material = nanotape_builder.get_material(nanotape_config)
        vacuum_config = VacuumConfiguration(
            size=vacuum_length,
            crystal=nanotape_material,
            direction=AxisEnum.x,
        )
        return NanoribbonConfiguration(
            stack_components=[nanotape_config, vacuum_config],
            direction=AxisEnum.x,
        ) 