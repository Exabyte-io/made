from typing import List, Union, Optional, Tuple

from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from mat3ra.esse.models.materials_category.pristine_structures.two_dimensional.slab import SlabConfigurationSchema

from mat3ra.made.material import Material
from ...stack.configuration import StackConfiguration
from ...vacuum.configuration import VacuumConfiguration
from .base_configurations import AtomicLayersUnique, AtomicLayersUniqueRepeatedConfiguration


class SlabConfiguration(SlabConfigurationSchema, StackConfiguration):
    type: str = "SlabConfiguration"
    stack_components: List[
        Union[AtomicLayersUnique, AtomicLayersUniqueRepeatedConfiguration, VacuumConfiguration]  # No Materials!
    ]
    direction: AxisEnum = AxisEnum.z

    @property
    def atomic_layers(self):
        return self.stack_components[0]

    @property
    def vacuum_configuration(self) -> VacuumConfiguration:
        return self.stack_components[1]

    @classmethod
    def from_parameters(
        cls,
        material_or_dict: Union[Material, dict],
        miller_indices: Tuple[int, int, int],
        number_of_layers: int,
        termination_formula: Optional[str] = None,
        vacuum: float = 10.0,
    ) -> "SlabConfiguration":
        from mat3ra.made.tools.analyze.lattice_planes import CrystalLatticePlanesMaterialAnalyzer
        from mat3ra.made.tools.build.slab.builders import AtomicLayersUniqueRepeatedBuilder
        from mat3ra.made.tools.build.slab.helpers import select_slab_termination

        if isinstance(material_or_dict, dict):
            material = Material.create(material_or_dict)
        else:
            material = material_or_dict

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
        return cls(
            stack_components=[atomic_layers_repeated_configuration, vacuum_configuration],
            direction=AxisEnum.z,
        ) 