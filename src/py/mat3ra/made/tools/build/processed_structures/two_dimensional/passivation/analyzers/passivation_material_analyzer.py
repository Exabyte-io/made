from typing import List

from ......analyze.slab import SlabMaterialAnalyzer
from ..enums import SurfaceTypesEnum


class PassivationMaterialAnalyzer(SlabMaterialAnalyzer):
    passivant: str = "H"
    bond_length: float = 1.0
    surface: SurfaceTypesEnum = SurfaceTypesEnum.TOP

    @property
    def passivant_coordinates(self) -> List[List[float]]:
        raise NotImplementedError(
            "This method should be implemented in subclasses to return passivant coordinates based on the surface type."
        )
