from mat3ra.esse.models.materials_category_components.entities.core.zero_dimensional.atom import AtomSchema

from ..defect_site.configuration import PointDefectSiteConfiguration
from ....vacuum.configuration import VacuumConfiguration


class AtomAtCoordinateConfiguration(VacuumConfiguration, PointDefectSiteConfiguration):
    element: AtomSchema
