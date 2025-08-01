from typing import List, Union, Optional, Tuple

from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from mat3ra.esse.models.materials_category.pristine_structures.two_dimensional.slab import SlabConfigurationSchema

from mat3ra.made.material import Material
from ... import MaterialWithBuildMetadata
from ...stack.configuration import StackConfiguration
from ...vacuum.configuration import VacuumConfiguration
from .base_configurations import AtomicLayersUniqueRepeatedConfiguration

from mat3ra.made.tools.analyze.lattice_planes import CrystalLatticePlanesMaterialAnalyzer
from mat3ra.made.tools.build.slab.termination_utils import select_slab_termination


class SlabConfiguration(SlabConfigurationSchema, StackConfiguration):
    type: str = "SlabConfiguration"
    stack_components: List[Union[AtomicLayersUniqueRepeatedConfiguration, VacuumConfiguration]]  # No Materials!
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

    def set_vacuum(self, vacuum: float) -> None:
        vacuum_configuration = self.vacuum_configuration
        vacuum_configuration.size = vacuum
        self.stack_components[1] = vacuum_configuration

    @classmethod
    def from_parameters(
        cls,
        material_or_dict: Union[Material, dict],
        miller_indices: Tuple[int, int, int],
        number_of_layers: int,
        termination_formula: Optional[str] = None,
        vacuum: float = 10.0,
        use_conventional_cell: bool = True,
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
            material = MaterialWithBuildMetadata.create(material_or_dict)
        else:
            material = material_or_dict

        crystal_lattice_planes_analyzer = CrystalLatticePlanesMaterialAnalyzer(
            material=material, miller_indices=miller_indices
        )
        terminations = crystal_lattice_planes_analyzer.terminations
        termination = select_slab_termination(terminations, termination_formula)

        if use_conventional_cell:
            material = crystal_lattice_planes_analyzer.material_with_conventional_lattice
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
