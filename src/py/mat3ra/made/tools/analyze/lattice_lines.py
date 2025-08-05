from typing import Tuple

from mat3ra.made.material import Material

from .lattice_planes import CrystalLatticePlanesMaterialAnalyzer


class CrystalLatticeLinesMaterialAnalyzer(CrystalLatticePlanesMaterialAnalyzer):
    """
    Analyzer for crystal lattice lines, used for nanoribbon creation.

    This analyzer treats the (u,v) Miller indices as a surface with (u,v,1) for pymatgen,
    allowing us to get terminations and shifts for 1D line structures.
    """

    miller_indices_2d: Tuple[int, int]

    def __init__(self, material: Material, miller_indices_2d: Tuple[int, int], **kwargs):
        miller_indices = (miller_indices_2d[0], miller_indices_2d[1], 1)
        super().__init__(
            material=material, miller_indices=miller_indices, miller_indices_2d=miller_indices_2d, **kwargs
        )

    @property
    def all_lines_as_pymatgen_slabs_with_vacuum(self):
        return self.all_planes_as_pymatgen_slabs_with_vacuum

    @property
    def all_lines_as_pymatgen_slabs_without_vacuum(self):
        return self.all_planes_as_pymatgen_slabs_without_vacuum
