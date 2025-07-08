from typing import Union

from mat3ra.code.entity import InMemoryEntityPydantic

from mat3ra.made.tools.analyze import BaseMaterialAnalyzer
from mat3ra.made.tools.build.metadata import MaterialMetadata
from mat3ra.made.tools.build.slab.builders import SlabWithAdditionalLayersBuilder
from mat3ra.made.tools.build.slab.configurations import SlabConfiguration, SlabWithAdditionalLayersConfiguration


class MatchedCellsSlabConfigurationHolder(InMemoryEntityPydantic):
    """
    Holds the original slab configuration and the slab configuration with additional layers.
    Used to match the c lattice parameter of the original slab with the slab with additional layers.
    """

    slab_with_additional_layers: SlabWithAdditionalLayersConfiguration
    slab_with_adjusted_vacuum: SlabConfiguration


class SlabMaterialAnalyzer(BaseMaterialAnalyzer):

    def get_slab_configuration(self) -> SlabConfiguration:
        metadata = MaterialMetadata(**self.material.metadata)
        slab_build_configuration_dict = metadata.build.configuration
        if slab_build_configuration_dict.get("type") != "SlabConfiguration":
            raise ValueError("Material is not a slab.")

        return SlabConfiguration(**slab_build_configuration_dict)

    def get_slab_with_additional_layers_configuration_holder(
        self,
        additional_layers: Union[int, float] = 1,
        vacuum_thickness: float = 5.0,
    ) -> MatchedCellsSlabConfigurationHolder:

        slab_configuration = self.get_slab_configuration()
        atomic_layers = slab_configuration.atomic_layers
        original_layers = atomic_layers.number_of_repetitions
        crystal = atomic_layers.crystal
        miller_indices = atomic_layers.miller_indices
        termination_top = atomic_layers.termination_top

        if isinstance(termination_top, dict):
            termination_formula = termination_top.get("chemical_elements")
        else:
            termination_formula = termination_top.formula

        configuration_with_added_layers = SlabWithAdditionalLayersConfiguration.from_parameters(
            material_or_dict=crystal,
            miller_indices=miller_indices,
            number_of_layers=original_layers,
            number_of_additional_layers=additional_layers,
            termination_formula=termination_formula,
            vacuum=vacuum_thickness,
        )

        slab_with_additional_layers = SlabWithAdditionalLayersBuilder().get_material(configuration_with_added_layers)
        # Calculate the vacuum thickness needed for the original slab to match the height of the slab with additional layers
        delta_c = slab_with_additional_layers.lattice.c - self.material.lattice.c

        vacuum_to_match_height = (
            SlabConfiguration(
                **MaterialMetadata(**self.material.metadata).build.configuration
            ).vacuum_configuration.size
            + delta_c
        )

        configuration_original_with_adjusted_vacuum = SlabConfiguration.from_parameters(
            material_or_dict=crystal,
            miller_indices=miller_indices,
            number_of_layers=original_layers,
            termination_formula=termination_formula,
            vacuum=vacuum_to_match_height,
        )
        return MatchedCellsSlabConfigurationHolder(
            slab_with_additional_layers=configuration_with_added_layers,
            slab_with_adjusted_vacuum=configuration_original_with_adjusted_vacuum,
        )
