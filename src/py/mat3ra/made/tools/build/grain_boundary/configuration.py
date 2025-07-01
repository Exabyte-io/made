from typing import Optional, List, Tuple

from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum

from mat3ra.made.material import Material
from ..interface.configuration import InterfaceConfiguration
from ..slab.configurations import SlabConfiguration
from ..vacuum.configuration import VacuumConfiguration


class GrainBoundaryConfiguration(InterfaceConfiguration):
    """
    Configuration for creating a grain boundary between two phases.
    """

    direction: AxisEnum = AxisEnum.x

    @property
    def phase_1_configuration(self) -> SlabConfiguration:
        return self.stack_components[0]

    @property
    def phase_2_configuration(self) -> SlabConfiguration:
        return self.stack_components[1]

    @property
    def translation_vector(self) -> List[float]:
        return [0.0, self.xy_shift[0], self.xy_shift[1]]

    @property
    def gap(self) -> float:
        return self.vacuum_configuration.size if len(self.stack_components) > 2 else 0.0

    @classmethod
    def from_parameters(
        cls,
        phase_1_material: Material,
        phase_2_material: Optional[Material] = None,
        phase_1_miller_indices: Tuple[int, int, int] = (0, 0, 1),
        phase_2_miller_indices: Tuple[int, int, int] = (0, 0, 1),
        phase_1_thickness: int = 1,
        phase_2_thickness: int = 1,
        xy_shift: Optional[List[float]] = None,
        gap: Optional[float] = None,
    ) -> "GrainBoundaryConfiguration":
        if phase_2_material is None:
            phase_2_material = phase_1_material
        if xy_shift is None:
            xy_shift = [0.0, 0.0]

        stack_components = [
            SlabConfiguration.from_material(
                material=phase_1_material,
                miller_indices=phase_1_miller_indices,
                number_of_layers=phase_1_thickness,
            ),
            SlabConfiguration.from_material(
                material=phase_2_material,
                miller_indices=phase_2_miller_indices,
                number_of_layers=phase_2_thickness,
                xy_shift=xy_shift,
            ),
        ]

        if gap is not None and gap > 0:
            stack_components.append(
                VacuumConfiguration(size=gap, crystal=phase_2_material.crystal, direction=AxisEnum.z)
            )

        return cls(stack_components=stack_components, direction=AxisEnum.x, xy_shift=xy_shift)


class GrainBoundaryWithVacuumConfiguration(SlabConfiguration, GrainBoundaryConfiguration):
    """
    Configuration for creating a grain boundary with vacuum (slab).

    This configuration is used to create a grain boundary slab with vacuum
    on both sides, similar to how SlabConfiguration works.
    """

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
        if xy_supercell_matrix is None:
            xy_supercell_matrix = [[1, 0], [0, 1]]

        return cls(
            stack_components=[
                SlabConfiguration.from_material(
                    material=grain_boundary_material,
                    miller_indices=miller_indices,
                    number_of_layers=number_of_layers,
                    vacuum=0.0,
                    termination_formula=termination_formula,
                    xy_supercell_matrix=xy_supercell_matrix,
                ),
                VacuumConfiguration(size=vacuum, crystal=grain_boundary_material.crystal, direction=AxisEnum.z),
            ],
            direction=AxisEnum.z,
        )
