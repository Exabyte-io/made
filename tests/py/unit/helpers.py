from typing import Optional, Tuple

from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.lattice_planes import CrystalLatticePlanesMaterialAnalyzer
from mat3ra.made.tools.build.slab.builders import AtomicLayersUniqueRepeatedBuilder
from mat3ra.made.tools.build.slab.configuration import (
    AtomicLayersUniqueRepeatedConfiguration,
    SlabConfiguration,
    VacuumConfiguration,
)
from mat3ra.made.tools.build.slab.helpers import select_slab_termination


def get_slab_configuration(
    material_dict: dict,
    miller_indices: Tuple[int, int, int],
    number_of_layers: int,
    termination_formula: Optional[str] = None,
    vacuum: float = 10.0,
) -> SlabConfiguration:
    """
    Helper function to create a SlabConfiguration.
    """
    material = Material.create(material_dict)
    crystal_lattice_planes_analyzer = CrystalLatticePlanesMaterialAnalyzer(
        material=material, miller_indices=miller_indices
    )
    terminations = crystal_lattice_planes_analyzer.terminations
    termination = select_slab_termination(terminations, termination_formula)

    atomic_layers_repeated_configuration = AtomicLayersUniqueRepeatedConfiguration(
        crystal=material,
        miller_indices=miller_indices,
        termination_top=termination,
        number_of_repetitions=number_of_layers,
    )
    atomic_layers_repeated_orthogonal_c = AtomicLayersUniqueRepeatedBuilder().get_material(
        atomic_layers_repeated_configuration
    )
    vacuum_configuration = VacuumConfiguration(
        size=vacuum, crystal=atomic_layers_repeated_orthogonal_c, direction=AxisEnum.z
    )
    return SlabConfiguration(
        stack_components=[atomic_layers_repeated_configuration, vacuum_configuration],
        direction=AxisEnum.z,
    )
