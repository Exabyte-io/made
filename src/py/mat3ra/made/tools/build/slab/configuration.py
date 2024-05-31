from typing import List, Tuple
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer as PymatgenSpacegroupAnalyzer

from mat3ra.made.material import Material
from ...convert import from_pymatgen, to_pymatgen, PymatgenStructure
from .builders import SlabBuildParameters
from .builders import SlabBuilder, SlabSelectorParameters


class BaseSlabConfiguration(object):
    @property
    def bulk(self) -> Material:
        raise NotImplementedError

    @property
    def miller_indices(self) -> Tuple[int, int, int]:
        raise NotImplementedError


class SlabConfiguration(BaseSlabConfiguration):
    """Class to generate slabs and manage different terminations."""

    thickness: int = 1
    vacuum: float = 0.5

    def __init__(
        self,
        bulk: Material = Material(Material.default_config),
        miller_indices: Tuple[int, int, int] = (0, 0, 1),
        thickness: int = 1,
        vacuum: float = 0.5,
        xy_supercell_matrix: List[List[int]] = [[1, 0], [0, 1]],
        use_conventional_cell: bool = True,
        use_orthogonal_z: bool = False,
    ):
        super().__init__()
        self.__bulk: PymatgenStructure = (
            PymatgenSpacegroupAnalyzer(to_pymatgen(bulk)).get_conventional_standard_structure()
            if use_conventional_cell
            else to_pymatgen(bulk)
        )
        self.__miller_indices = miller_indices
        self.thickness = thickness
        self.vacuum = vacuum
        self.__builder = SlabBuilder(
            build_parameters=SlabBuildParameters(
                use_orthogonal_z=use_orthogonal_z, xy_supercell_matrix=xy_supercell_matrix
            )
        )

    @property
    def bulk(self):
        return from_pymatgen(self.__bulk)

    @property
    def miller_indices(self):
        return self.__miller_indices

    def get_slab(self, termination) -> Material:
        return self.__builder.get_material(self, selector_parameters=SlabSelectorParameters(termination=termination))

    @property
    def terminations(self):
        return self.__builder.terminations(self)
