from math import ceil
from typing import Any

from pydantic import BaseModel

from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.slab import SlabMaterialAnalyzer
from mat3ra.made.tools.build.defect.slab.configuration import SlabDefectConfiguration, SlabStackConfiguration
from mat3ra.made.tools.build.merge import MergeBuilder
from mat3ra.made.tools.build.metadata import MaterialMetadata
from mat3ra.made.tools.build.slab.builders import SlabBuilder
from mat3ra.made.tools.build.slab.configurations import SlabConfiguration
from mat3ra.made.tools.build.slab.helpers import create_slab
from mat3ra.made.tools.build.stack.builders import StackNComponentsBuilder
from mat3ra.made.tools.modify import filter_by_box


def get_slab_build_configuration(metadata: dict) -> SlabConfiguration:

    material_metadata = MaterialMetadata(**metadata)
    slab_build_configuration_dict = material_metadata.current_build.configuration

    if slab_build_configuration_dict.get("type") != "SlabConfiguration":
        raise ValueError("Material is not a slab.")

    return SlabConfiguration(**slab_build_configuration_dict)


class SlabDefectBuilderParameters(BaseModel):
    auto_add_vacuum: bool = True
    vacuum: float = 5.0


class SlabDefectBuilder(MergeBuilder):
    _ConfigurationType = SlabDefectConfiguration

    def get_slab_with_additional_layers(
        self, slab: Material, number_of_additional_layers: int = 1, vacuum: float = 0.0
    ) -> Material:
        """
        Create a slab with additional layers by merging a slab material with additional layers and vacuum.
        Args:
            slab: The original slab material.
            number_of_additional_layers: Number of additional layers to add to the slab.
            vacuum: Thickness of the vacuum layer to be added.
        Returns:
            Material: The new slab material with additional layers and vacuum if needed.
        """
        # get_slab_build_configuration handles the supercell creation, etc.
        slab_configuration = get_slab_build_configuration(slab.metadata)
        total_number_of_layers = slab_configuration.number_of_layers + number_of_additional_layers
        return create_slab(
            crystal=slab,
            miller_indices=slab_configuration.miller_indices,
            termination=slab_configuration.termination,
            number_of_layers=total_number_of_layers,
            vacuum=vacuum,
            xy_supercell_matrix=slab_configuration.xy_supercell_matrix,
        )

    def get_slab_with_additional_layers_float(
        self, slab: Material, number_of_additional_layers: float = 0.5, vacuum: float = 0.0
    ) -> Material:
        """
        Create a slab with additional layers by merging a slab material with additional layers and vacuum.
        Args:
            slab: The original slab material.
            number_of_additional_layers: Number of additional layers to add to the slab.
            vacuum: Thickness of the vacuum layer to be added.
        Returns:
            Material: The new slab material with additional layers and vacuum if needed.
        """
        slab_analyzer = SlabMaterialAnalyzer(material=slab)
        ceiling_number_of_additional_layers = int(ceil(number_of_additional_layers))
        slab_with_int_additional_layers = self.get_slab_with_additional_layers(
            slab=slab,
            number_of_additional_layers=ceiling_number_of_additional_layers,
            vacuum=vacuum,
        )
        vacuum_thickness_in_layers = slab_analyzer.vacuum_thickness_in_layers
        number_of_layers = slab_analyzer.number_of_layers
        max_z_crystal_coordinate = (number_of_additional_layers + number_of_layers + vacuum_thickness_in_layers) / (
            ceiling_number_of_additional_layers + number_of_layers + vacuum_thickness_in_layers
        )
        return filter_by_box(
            slab_with_int_additional_layers,
            max_coordinate=[1, 1, max_z_crystal_coordinate],
            reset_ids=True,
        )


class SlabStackBuilder(StackNComponentsBuilder):
    def _configuration_to_material(self, configuration_or_material: Any) -> Material:
        if isinstance(configuration_or_material, SlabConfiguration):
            return SlabBuilder().get_material(configuration_or_material)
        else:
            return super()._configuration_to_material(configuration_or_material)

    def _generate(self, configuration: SlabStackConfiguration) -> Material:
        return super()._generate(configuration)
