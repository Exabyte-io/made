from typing import Optional, Any

from mat3ra.made.tools.build_components import MaterialWithBuildMetadata
from ..base.builder import PointDefectBuilder
from .configuration import VacancyDefectConfiguration


class VacancyDefectBuilder(PointDefectBuilder):
    _ConfigurationType = VacancyDefectConfiguration

    def _post_process(
        self,
        material: MaterialWithBuildMetadata,
        post_process_parameters: Optional[Any] = None,
        configuration: Optional[TypeConfiguration] = None,
    ) -> MaterialWithBuildMetadata:
        material = super()._post_process(material, post_process_parameters, configuration)
        material.basis.remove_atoms_by_elements("Vac")
        return material
