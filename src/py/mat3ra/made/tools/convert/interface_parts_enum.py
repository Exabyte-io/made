from enum import Enum


class InterfacePartsEnum(int, Enum):
    SUBSTRATE = 0
    FILM = 1


INTERFACE_LABELS_MAP = {"substrate": InterfacePartsEnum.SUBSTRATE.value, "film": InterfacePartsEnum.FILM.value}
