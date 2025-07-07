from typing import Union, Tuple

from mat3ra.made.tools.analyze import BaseMaterialAnalyzer
from mat3ra.made.tools.analyze.other import get_atomic_coordinates_extremum
from mat3ra.made.tools.build.metadata import MaterialMetadata
from mat3ra.made.tools.build.slab.configurations import SlabConfiguration
from mat3ra.made.tools.build.slab.configurations import SlabWithAdditionalLayersConfiguration


class SlabMaterialAnalyzer(BaseMaterialAnalyzer):
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
        metadata = MaterialMetadata(**self.material.metadata)
        slab_build_configuration_dict = metadata.build.configuration
        if slab_build_configuration_dict.get("type") != "SlabConfiguration":
            raise ValueError("Material is not a slab.")

        adjusted_vacuum = self._calculate_required_vacuum(additional_layers, vacuum_thickness)

        stack_components = slab_build_configuration_dict.get("stack_components", [])
        if not stack_components or len(stack_components) < 2:
            raise ValueError("Invalid slab configuration: missing stack components")

        atomic_layers_config = stack_components[0]
        original_layers = atomic_layers_config.get("number_of_repetitions", 1)
        miller_indices = atomic_layers_config.get("miller_indices", [0, 0, 1])
        termination_formula = atomic_layers_config.get("termination_top", "Si")
        crystal = atomic_layers_config.get("crystal", self.material)

        configuration_original_with_adjusted_vacuum = SlabConfiguration.from_parameters(
            material_or_dict=crystal,
            miller_indices=miller_indices,
            number_of_layers=original_layers,
            termination_formula=termination_formula,
            vacuum=adjusted_vacuum,
        )

        configuration_with_added_layers = SlabWithAdditionalLayersConfiguration.from_parameters(
            material_or_dict=crystal,
            miller_indices=miller_indices,
            number_of_layers=original_layers,
            additional_layers=additional_layers,
            termination_formula=termination_formula,
            vacuum=adjusted_vacuum,
        )

        return (configuration_with_added_layers, configuration_original_with_adjusted_vacuum)

    def _calculate_required_vacuum(self, additional_layers: Union[int, float], vacuum_thickness: float) -> float:
        layer_height = self._get_layer_height_cartesian()
        additional_height = additional_layers * layer_height
        current_max_z_crystal = get_atomic_coordinates_extremum(
            self.material, "max", "z", use_cartesian_coordinates=False
        )
        new_max_z_crystal = current_max_z_crystal + additional_height
        if new_max_z_crystal > 1.0:
            required_vacuum = (new_max_z_crystal - 1.0) * self.material.lattice.c
            return max(required_vacuum, vacuum_thickness)

        return vacuum_thickness

    def _get_layer_height_cartesian(self) -> float:
        metadata = MaterialMetadata(**self.material.metadata)
        build_config = metadata.build.configuration
        original_layers = build_config.get("number_of_layers", 1)
        current_max_z_crystal = get_atomic_coordinates_extremum(
            self.material, "max", "z", use_cartesian_coordinates=False
        )
        layer_height = current_max_z_crystal / original_layers

        return layer_height
