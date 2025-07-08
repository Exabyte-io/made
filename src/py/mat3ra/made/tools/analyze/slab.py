from typing import Tuple, Union

from mat3ra.made.tools.analyze import BaseMaterialAnalyzer
from mat3ra.made.tools.analyze.other import get_atomic_coordinates_extremum
from mat3ra.made.tools.build.metadata import MaterialMetadata
from mat3ra.made.tools.build.slab.builders import SlabBuilder, SlabWithAdditionalLayersBuilder
from mat3ra.made.tools.build.slab.configurations import SlabConfiguration, SlabWithAdditionalLayersConfiguration


class SlabMaterialAnalyzer(BaseMaterialAnalyzer):

    def get_slab_configuration(self) -> SlabConfiguration:
        metadata = MaterialMetadata(**self.material.metadata)
        slab_build_configuration_dict = metadata.build.configuration
        if slab_build_configuration_dict.get("type") != "SlabConfiguration":
            raise ValueError("Material is not a slab.")

        return SlabConfiguration(**slab_build_configuration_dict)

    def get_slab_with_additional_layers_configurations(
        self,
        additional_layers: Union[int, float] = 1,
        vacuum_thickness: float = 5.0,
    ) -> Tuple[SlabWithAdditionalLayersConfiguration, SlabConfiguration]:
        """
        Analyze a slab material and create two configurations:
        - with original number of layers and adjusted vacuum to match the new height
        - with additional layers and the vacuum thickness on top.
        Returns:
            (SlabWithAdditionalLayersConfiguration, SlabConfiguration)
        """
        slab_configuration = self.get_slab_configuration()
        atomic_layers = slab_configuration.atomic_layers
        original_layers = atomic_layers.number_of_repetitions
        miller_indices = atomic_layers.miller_indices
        termination_top = atomic_layers.termination_top
        crystal = atomic_layers.crystal

        if isinstance(termination_top, dict):
            termination_formula = termination_top.get("chemical_elements")
        else:
            termination_formula = termination_top.formula

        # Calculate the vacuum thickness needed for the original slab to match the height of the slab with additional layers
        vacuum_to_match_height = self._vacuum_thickness_for_original_slab(
            additional_layers=additional_layers,
            vacuum_thickness=vacuum_thickness,
        )

        configuration_with_added_layers = SlabWithAdditionalLayersConfiguration.from_parameters(
            material_or_dict=crystal,
            miller_indices=miller_indices,
            number_of_layers=original_layers,
            number_of_additional_layers=additional_layers,
            termination_formula=termination_formula,
            vacuum=vacuum_thickness,
        )
        configuration_original_with_adjusted_vacuum = SlabConfiguration.from_parameters(
            material_or_dict=crystal,
            miller_indices=miller_indices,
            number_of_layers=original_layers,
            termination_formula=termination_formula,
            vacuum=vacuum_to_match_height,
        )

        return (configuration_with_added_layers, configuration_original_with_adjusted_vacuum)

    def _vacuum_thickness_for_original_slab(
        self, additional_layers: Union[int, float] = 1, vacuum_thickness: float = 5.0
    ) -> float:
        """
        Calculate the vacuum thickness needed for the original slab to match the height of the slab with additional layers.

        The approach is to calculate the height that will be added by additional layers.
        Each additional layer adds exactly the lattice parameter 'a' to the c vector.
        """
        slab_configuration = self.get_slab_configuration()

        slab = SlabBuilder().get_material(slab_configuration)
        lattice_parameter_a = slab.lattice.vector_arrays[0][0]

        additional_height = lattice_parameter_a * additional_layers

        return additional_height + vacuum_thickness
