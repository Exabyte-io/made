from typing import List, Optional

import numpy as np
from mat3ra.code.array_with_ids import ArrayWithIds
from pydantic import BaseModel


class CrystalSite(BaseModel):
    # element: str
    # coordinate: List[float]
    nearest_neighbor_vectors: List[np.ndarray] = []
    # coordination_number: int = 0
    # see https://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-wp-list for an example
    wyckoff_letter: Optional[str] = None

    class Config:
        arbitrary_types_allowed = True

    @property
    def coordination_number(self):
        return len(self.nearest_neighbor_vectors)


class CrystalSiteList(ArrayWithIds):
    values: List[CrystalSite]
