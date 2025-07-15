from typing import List

from mat3ra.made.tools.analyze.slab import SlabMaterialAnalyzer


class NanoparticleMaterialAnalyzer(SlabMaterialAnalyzer):
    orientation_z = [0, 0, 1]
    vacuum_padding = 10.0

    @property
    def supercell_to_fit_nanpoarticle(self) -> List[List[int]]:
        """
        Calculate the supercell size needed to fit the nanoparticle and vacuum padding.
        """
