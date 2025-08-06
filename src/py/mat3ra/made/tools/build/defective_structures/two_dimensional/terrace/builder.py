from typing import Type

from .build_parameters import TerraceBuildParameters
from .configuration import TerraceDefectConfiguration
from .....build_components import MaterialWithBuildMetadata
from .....build_components.entities.reusable.two_dimensional.slab_stack.builder import SlabStackBuilder
from .....modify import translate_to_z_level
from .....operations.core.unary import edit_cell, rotate


class TerraceDefectBuilder(SlabStackBuilder):
    _BuildParametersType: Type[TerraceBuildParameters] = TerraceBuildParameters

    def get_name_suffix(self, configuration: TerraceDefectConfiguration) -> str:
        return f"Terrace {configuration.cut_direction}"

    def _generate(self, configuration: TerraceDefectConfiguration) -> MaterialWithBuildMetadata:
        material = super()._generate(configuration)
        if self.build_parameters.rotate_to_match_pbc:
            axis = self.build_parameters.axis
            angle = self.build_parameters.angle
            new_lattice_vectors = self.build_parameters.new_lattice_vectors
            translated_material = translate_to_z_level(material, "center")
            adjusted_material = edit_cell(translated_material, new_lattice_vectors)
            material = rotate(adjusted_material, axis, angle)

        return material
