from typing import List, Union, Optional, Tuple

from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from mat3ra.esse.models.materials_category.pristine_structures.two_dimensional.slab import SlabConfigurationSchema

from mat3ra.made.material import Material
from ...stack.configuration import StackConfiguration
from ...vacuum.configuration import VacuumConfiguration
from .base_configurations import AtomicLayersUniqueConfiguration, AtomicLayersUniqueRepeatedConfiguration

from mat3ra.made.tools.analyze.lattice_planes import CrystalLatticePlanesMaterialAnalyzer
from mat3ra.made.tools.build.slab.termination_utils import select_slab_termination


class SlabConfiguration(SlabConfigurationSchema, StackConfiguration):
    type: str = "SlabConfiguration"
    stack_components: List[
        Union[
            AtomicLayersUniqueConfiguration, AtomicLayersUniqueRepeatedConfiguration, VacuumConfiguration
        ]  # No Materials!
    ]
    direction: AxisEnum = AxisEnum.z

    @property
    def number_of_layers(self):
        return self.atomic_layers.number_of_repetitions

    @property
    def vacuum(self):
        return self.vacuum_configuration.size

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
        """
        Creates a SlabConfiguration from the given parameters.
        Args:
            material_or_dict (Union[Material, dict]): Material or dictionary representation of the material.
            miller_indices (Tuple[int, int, int]): Miller indices for the slab surface.
            number_of_layers (int): Number of atomic layers in the slab, in the number of unit cells.
            termination_formula (Optional[str]): Formula of the termination to use for the slab (i.e. "SrTiO").
            vacuum (float): Size of the vacuum layer in Angstroms.

        Returns:
            SlabConfiguration: The created slab configuration.
        """
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

        vacuum_configuration = VacuumConfiguration(size=vacuum, crystal=None, direction=AxisEnum.z)

        return cls(
            stack_components=[atomic_layers_repeated_configuration, vacuum_configuration],
            direction=AxisEnum.z,
        )

    def to_parameters(self) -> dict:
        atomic_layers = self.atomic_layers
        return {
            "material_or_dict": atomic_layers.crystal,
            "miller_indices": atomic_layers.miller_indices,
            "number_of_layers": atomic_layers.number_of_repetitions,
            "termination_formula": atomic_layers.termination_top.formula if atomic_layers.termination_top else None,
            "vacuum": self.vacuum_configuration.size,
        }
