from mat3ra.esse.models.materials_category_components.entities.core.zero_dimensional.atom import (
    AtomSchema,
)

from ...two_dimensional.vacuum.configuration import VacuumConfiguration
from ....auxiliary.zero_dimensional.point_defect_site.configuration import PointDefectSiteConfiguration


class AtomAtCoordinateConfiguration(VacuumConfiguration, PointDefectSiteConfiguration):
    element: AtomSchema
