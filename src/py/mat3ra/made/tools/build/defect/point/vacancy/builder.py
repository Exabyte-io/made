from typing import Optional, Any

from mat3ra.made.tools.build import MaterialWithBuildMetadata, TConfiguration
from mat3ra.made.tools.build.defect.point.base.builder import PointDefectBuilder
from mat3ra.made.tools.build.defect.point.vacancy.configuration import VacancyDefectConfiguration


class VacancyDefectBuilder(PointDefectBuilder):
    _ConfigurationType = VacancyDefectConfiguration

    def _post_process(
        self,
        material: MaterialWithBuildMetadata,
        post_process_parameters: Optional[Any] = None,
        configuration: Optional[TConfiguration] = None,
    ) -> MaterialWithBuildMetadata:
        material = super()._post_process(material, post_process_parameters, configuration)
        material.basis.remove_atoms_by_elements("Vac")
        return material
