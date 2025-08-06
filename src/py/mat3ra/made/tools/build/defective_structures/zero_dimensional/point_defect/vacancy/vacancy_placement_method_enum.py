from enum import Enum

from .build.defective_structures.zero_dimensional.point_defect.atom_placement_method_enum import \
    AtomPlacementMethodEnum


class VacancyPlacementMethodEnum(Enum):
    CLOSEST_SITE = AtomPlacementMethodEnum.CLOSEST_SITE
