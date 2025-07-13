from mat3ra.made.material import Material
from ... import MaterialWithBuildMetadata
from ..slab.builders import SlabStackBuilder
from ..slab.configuration import SlabStackConfiguration
from .configuration import TerraceDefectConfiguration
from ...defect.terrace.parameters import TerraceBuildParameters
from mat3ra.made.tools.modify import rotate


class TerraceDefectBuilder(SlabStackBuilder):
    _BuildParametersType = TerraceBuildParameters

    def _generate(self, configuration: TerraceDefectConfiguration) -> Material:
        material = super()._generate(configuration)
        if self.build_parameters.rotate_to_match_pbc:
            axis = self.build_parameters.axis
            angle = self.build_parameters.angle
            material = rotate(material, axis, angle)

        return material

    def _update_material_name(
        self, material: MaterialWithBuildMetadata, configuration: SlabStackConfiguration
    ) -> MaterialWithBuildMetadata:
        new_material = super()._update_material_name(material, configuration)
        new_name = f"{new_material.name}, Terrace {configuration.cut_direction}"
        new_material.name = new_name
        return new_material
