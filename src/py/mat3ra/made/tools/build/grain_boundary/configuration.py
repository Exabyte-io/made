from typing import Optional, List, Tuple

from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.interface.utils.holders import MatchedSubstrateFilmConfigurationHolder
from .. import BaseConfiguration, BaseConfigurationPydantic
from ..slab.configurations import SlabConfiguration


class GrainBoundaryConfiguration(BaseConfigurationPydantic):
    """
    Configuration for creating a grain boundary between two phases.

    The grain boundary is created by stacking two phases in x/y direction
    with a relative shift perpendicular to the interface.
    """

    phase_1_configuration: MatchedSubstrateFilmConfigurationHolder
    phase_2_configuration: MatchedSubstrateFilmConfigurationHolder
    translation_vector: List[float] = [0.0, 0.0, 0.0]  # Relative shift between phases
    gap: float = 3.0  # Gap between phases in Angstroms

    @property
    def substrate_configuration(self) -> SlabConfiguration:
        """Get the substrate configuration (phase 1)."""
        return self.phase_1_configuration.substrate_configuration

    @property
    def film_configuration(self) -> SlabConfiguration:
        """Get the film configuration (phase 2)."""
        return self.phase_2_configuration.film_configuration


class GrainBoundaryWithVacuumConfiguration(BaseConfigurationPydantic):
    """
    Configuration for creating a grain boundary with vacuum (slab).

    This configuration is used to create a grain boundary slab with vacuum
    on both sides, similar to how SlabConfiguration works.
    """

    grain_boundary_material: Material
    miller_indices: Tuple[int, int, int]
    number_of_layers: int = 1
    vacuum: float = 10.0
    termination_formula: Optional[str] = None
    xy_supercell_matrix: List[List[int]] = [[1, 0], [0, 1]]

    @classmethod
    def from_grain_boundary_material(
        cls,
        grain_boundary_material: Material,
        miller_indices: Tuple[int, int, int],
        number_of_layers: int = 1,
        vacuum: float = 10.0,
        termination_formula: Optional[str] = None,
        xy_supercell_matrix: Optional[List[List[int]]] = None,
    ) -> "GrainBoundaryWithVacuumConfiguration":
        """Create configuration from a grain boundary material."""
        if xy_supercell_matrix is None:
            xy_supercell_matrix = [[1, 0], [0, 1]]

        return cls(
            grain_boundary_material=grain_boundary_material,
            miller_indices=miller_indices,
            number_of_layers=number_of_layers,
            vacuum=vacuum,
            termination_formula=termination_formula,
            xy_supercell_matrix=xy_supercell_matrix,
        )
