from enum import Enum


class InterfacePartsEnum(str, Enum):
    SUBSTRATE = 0
    FILM = 1


INTERFACE_LABELS_MAP = {"substrate": InterfacePartsEnum.SUBSTRATE, "film": InterfacePartsEnum.FILM}
