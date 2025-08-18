from typing import List

from mat3ra.code.array_with_ids import ArrayWithIds

from .crystal_site import CrystalSite


class CrystalSiteList(ArrayWithIds):
    values: List[CrystalSite]
