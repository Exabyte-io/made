from typing import Type

from .build_parameters import TerraceBuildParameters
from .configuration import TerraceDefectConfiguration
from .....build_components import MaterialWithBuildMetadata
from .....build_components.entities.reusable.two_dimensional.slab_stack.builder import SlabStackBuilder
from .....modify import translate_to_z_level
from .....operations.core.unary import edit_cell, rotate


class TerraceDefectBuilder(SlabStackBuilder):
    _ConfigurationType: Type[TerraceDefectConfiguration] = TerraceDefectConfiguration
    _BuildParametersType: Type[TerraceBuildParameters] = TerraceBuildParameters
    _DefaultBuildParameters: TerraceBuildParameters = TerraceBuildParameters()

    def get_name_suffix(self, configuration: TerraceDefectConfiguration) -> str:
        return f"Terrace {configuration.cut_direction}"

    def _generate(self, configuration: TerraceDefectConfiguration) -> MaterialWithBuildMetadata:
        material = super()._generate(configuration)
        build_params = self.build_parameters if self.build_parameters else self._DefaultBuildParameters
        if build_params.rotate_to_match_pbc:
            axis = build_params.axis
            angle = build_params.angle
            new_lattice_vectors = build_params.new_lattice_vectors
            translated_material = translate_to_z_level(material, "center")
            adjusted_material = edit_cell(translated_material, new_lattice_vectors)
            material = rotate(adjusted_material, axis, angle)

        return material
