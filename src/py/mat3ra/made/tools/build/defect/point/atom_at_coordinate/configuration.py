from mat3ra.esse.models.materials_category_components.entities.core.zero_dimensional.atom import AtomSchema

from mat3ra.made.tools.build.defect.point.defect_site.configuration import PointDefectSiteConfiguration
from mat3ra.made.tools.build.vacuum.configuration import VacuumConfiguration


class AtomAtCoordinateConfiguration(VacuumConfiguration, PointDefectSiteConfiguration):
    element: AtomSchema
