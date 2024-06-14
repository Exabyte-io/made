from typing import List


class UnitCell:
    ax: float
    ay: float
    az: float
    bx: float
    by: float
    bz: float
    cx: float
    cy: float
    cz: float
    units: str

    def __init__(self, vectors):
        self.ax, self.ay, self.az, self.bx, self.by, self.bz, self.cx, self.cy, self.cz, self.units = vectors

    @property
    def vectorA(self) -> List[float]:
        return [self.ax, self.ay, self.az]

    @property
    def vectorB(self) -> List[float]:
        return [self.bx, self.by, self.bz]

    @property
    def vectorC(self) -> List[float]:
        return [self.cx, self.cy, self.cz]

    @property
    def axes(self) -> List[List[float]]:
        return [self.vectorA, self.vectorB, self.vectorC]
