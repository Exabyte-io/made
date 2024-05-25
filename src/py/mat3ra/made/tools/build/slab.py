from pymatgen.core.structure import Structure
from pymatgen.core.surface import SlabGenerator as PymatgenSlabGenerator
from pymatgen.core.interface import label_termination
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from typing import List, Tuple
from .supercell import create_supercell
from ..convert import decorator_convert_material_args_kwargs_to_structure, from_pymatgen
from ...material import Material


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
    xy_supercell_matrix: List[List[int]] = [[1, 0], [0, 1]]
    use_orthogonal_z: bool = False
    use_conventional_cell: bool = True

    @decorator_convert_material_args_kwargs_to_structure
    def __init__(
        self,
        bulk: Structure,
        miller_indices: Tuple[int, int, int] = (0, 0, 1),
        thickness: int = 1,
        vacuum: float = 0.5,
        xy_supercell_matrix: List[List[int]] = [[1, 0], [0, 1]],
        use_conventional_cell: bool = True,
        use_orthogonal_z: bool = False,
    ):
        super().__init__()
        self.__bulk = (
            SpacegroupAnalyzer(bulk).get_conventional_standard_structure() if self.use_conventional_cell else bulk
        )
        self.__miller_indices = miller_indices
        self.thickness = thickness
        self.vacuum = vacuum
        self.xy_supercell_matrix = xy_supercell_matrix
        self.use_conventional_cell = use_conventional_cell
        self.use_orthogonal_z = use_orthogonal_z

    @property
    def bulk(self):
        return from_pymatgen(self.__bulk)

    @property
    def miller_indices(self):
        return self.__miller_indices

    @property
    def generator(self):
        return PymatgenSlabGenerator(
            initial_structure=self.__bulk,
            miller_index=self.miller_indices,
            min_slab_size=self.thickness,
            min_vacuum_size=self.vacuum,
            in_unit_planes=True,
            reorient_lattice=True,
        )

    @property
    def __slabs_with_unique_terminations(self):
        return [slab.get_orthogonal_c_slab() if self.use_orthogonal_z else slab for slab in self.generator.get_slabs()]

    @property
    def terminations(self):
        return list(set(label_termination(slab) for slab in self.__slabs_with_unique_terminations))

    def get_material(self, termination: str) -> Material:
        for slab in self.__slabs_with_unique_terminations:
            if label_termination(slab) == termination:
                return create_supercell(
                    Material(from_pymatgen(slab)),
                    self.xy_supercell_matrix,
                )
        raise ValueError(f"Termination {termination} not found in slabs.")
