from mat3ra.esse.models.materials_category_components.entities.core.zero_dimensional.atom import AtomSchema

from ....auxiliary.zero_dimensional.point_defect_site.configuration import PointDefectSiteConfiguration
from ...two_dimensional.vacuum.configuration import VacuumConfiguration


class AtomAtCoordinateConfiguration(VacuumConfiguration, PointDefectSiteConfiguration):
    element: AtomSchema
