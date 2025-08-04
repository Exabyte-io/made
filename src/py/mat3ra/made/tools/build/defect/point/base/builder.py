from typing import Type, Dict, Union

from mat3ra.made.material import Material
from mat3ra.made.tools.build import MaterialWithBuildMetadata, TConfiguration
from .configuration import PointDefectConfiguration
from mat3ra.made.tools.build.defect.point.defect_site.builder import PointDefectSiteBuilder
from mat3ra.made.tools.build.defect.point.defect_site.configuration import PointDefectSiteConfiguration
from mat3ra.made.tools.build.defect.point.intersitital.configuration import (
    InterstitialDefectConfiguration,
)
from mat3ra.made.tools.build.defect.point.substitutional.configuration import SubstitutionalDefectConfiguration
from mat3ra.made.tools.build.defect.point.vacancy.configuration import VacancyDefectConfiguration
from mat3ra.made.tools.build.merge import MergeBuilder


class PointDefectBuilder(MergeBuilder):
    _ConfigurationType: Type[PointDefectConfiguration] = PointDefectConfiguration

    @property
    def merge_component_types_conversion_map(self) -> Dict[Type, Type]:
        return {
            VacancyDefectConfiguration: PointDefectSiteBuilder,
            SubstitutionalDefectConfiguration: PointDefectSiteBuilder,
            InterstitialDefectConfiguration: PointDefectSiteBuilder,
            PointDefectSiteConfiguration: PointDefectSiteBuilder,
        }

    def _update_material_name(
        self, material: Union[Material, MaterialWithBuildMetadata], configuration: TConfiguration
    ) -> MaterialWithBuildMetadata:
        host_material = None
        for component in configuration.merge_components:
            if isinstance(component, Material):
                host_material = component
                break

        if host_material:
            defect_type = configuration.__class__.__name__.lower()
            defect_type = defect_type.replace("defectconfiguration", "").replace("configuration", "")
            material.name = f"{host_material.name} with {defect_type} defect"

        return material
