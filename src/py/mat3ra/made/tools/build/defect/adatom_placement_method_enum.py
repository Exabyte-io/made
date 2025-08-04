from enum import Enum

from .atom_placement_method_enum import AtomPlacementMethodEnum


class AdatomPlacementMethodEnum(Enum):
    EXACT_COORDINATE = AtomPlacementMethodEnum.EXACT_COORDINATE
    EQUIDISTANT = AtomPlacementMethodEnum.EQUIDISTANT
    NEW_CRYSTAL_SITE = AtomPlacementMethodEnum.NEW_CRYSTAL_SITE