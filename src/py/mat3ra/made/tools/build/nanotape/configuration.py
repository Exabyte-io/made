from typing import List, Union, Optional, Tuple

from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum

from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.lattice_lines import CrystalLatticeLinesAnalyzer
from ..lattice_lines.builders import CrystalLatticeLinesRepeatedBuilder
from ..lattice_lines.configuration import (
    CrystalLatticeLinesUniqueRepeatedConfiguration,
    get_miller_indices_from_edge_type,
    EdgeTypes,
)
from ..slab.entities import Termination
from ..stack.configuration import StackConfiguration
from ..vacuum.configuration import VacuumConfiguration


class NanoTapeConfiguration(StackConfiguration):
    """
    Configuration for building a nanotape from crystal lattice lines.
    NanoTape = [CLLUR, vacuum] stacked on Y direction.

    Args:
        stack_components: List of configuration objects for nanotape components.
        direction: Direction along which to stack components.
    """

    type: str = "NanoTapeConfiguration"
    stack_components: List[Union["CrystalLatticeLinesUniqueRepeatedConfiguration", VacuumConfiguration]]
    direction: AxisEnum = AxisEnum.y
    use_rectangular_lattice: bool = True

    @property
    def lattice_lines(self):
        """Get the lattice lines configuration component."""
        return self.stack_components[0]

    @property
    def vacuum_configuration(self) -> VacuumConfiguration:
        """Get the vacuum configuration component."""
        return self.stack_components[1]

    @classmethod
    def from_parameters(
        cls,
        material: Material,
        miller_indices_uv: Optional[Tuple[int, int]] = None,
        edge_type: Optional[EdgeTypes] = EdgeTypes.zigzag,
        width: int = 2,
        length: int = 2,
        vacuum_width: float = 10.0,
        termination: Optional[Termination] = None,
    ) -> "NanoTapeConfiguration":
        """
        Create a NanoTapeConfiguration from parameters.

        Args:
            material: The monolayer material to create the nanotape from.
            miller_indices_uv: The (u,v) Miller indices for the nanotape direction.
            edge_type: Edge type string ("zigzag"/"armchair"). Optional if miller_indices_uv is provided.
            width: The width of the nanotape in number of unit cells.
            length: The length of the nanotape in number of unit cells.
            vacuum_width: The width of the vacuum region in Angstroms (cartesian).
            termination: The termination to use for the nanotape. If None, uses default termination.

        Returns:
            NanoTapeConfiguration: The nanotape configuration.
        """
        if miller_indices_uv is None and edge_type is None:
            raise ValueError("Either miller_indices_uv or edge_type must be provided")
        if miller_indices_uv is None and edge_type is not None:
            miller_indices_uv = get_miller_indices_from_edge_type(edge_type)

        lattice_lines_analyzer = CrystalLatticeLinesAnalyzer(material=material, miller_indices_uv=miller_indices_uv)
        if termination is None:
            termination = lattice_lines_analyzer.default_termination

        lattice_lines_config = CrystalLatticeLinesUniqueRepeatedConfiguration(
            crystal=material,
            miller_indices_uv=miller_indices_uv,
            termination_top=termination,
            number_of_repetitions_width=width,
            number_of_repetitions_length=length,
        )

        lattice_lines_builder = CrystalLatticeLinesRepeatedBuilder()
        lattice_lines_material = lattice_lines_builder.get_material(lattice_lines_config)

        vacuum_config = VacuumConfiguration(
            size=vacuum_width,
            crystal=lattice_lines_material,
            direction=AxisEnum.y,
        )

        return cls(
            stack_components=[lattice_lines_config, vacuum_config],
            direction=AxisEnum.y,
        )
