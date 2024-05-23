from pymatgen.core.structure import Structure
from pymatgen.core.surface import SlabGenerator as PymatgenSlabGenerator
from pymatgen.core.interface import label_termination
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from typing import Any, List, Tuple, Dict, Optional
import json
from pydantic import BaseModel
from .supercell import create_supercell
from ..convert import decorator_convert_material_args_kwargs_to_structure, from_pymatgen
from ...material import Material


# @decorator_convert_material_args_kwargs_to_structure
class BaseSlabConfiguration(BaseModel):
    bulk: Structure
    miller_indices: Tuple[int, int, int] = (0, 0, 1)

    @property
    def bulk(self):
        return self.__bulk

    @property
    def miller_indices(self):
        return self.__miller_indices


class SlabConfiguration(BaseSlabConfiguration):
    """Class to generate slabs and manage different terminations."""

    @decorator_convert_material_args_kwargs_to_structure
    def __init__(
        self,
        bulk: Structure,
        miller_indices: Tuple[int, int, int] = (0, 0, 1),
        thickness: int = 1,
        vacuum: float = 0.5,
        xy_supercell_matrix: List[List[int]] = [[1, 0], [0, 1]],
        use_conventional_cell: bool = True,
    ):
        super().__init__(bulk=bulk, miller_indices=miller_indices)
        self.__bulk = bulk
        self.__miller_indices = miller_indices
        self.thickness = thickness
        self.vacuum = vacuum
        self.xy_supercell_matrix = xy_supercell_matrix
        self.use_conventional_cell = use_conventional_cell

    @property
    def bulk(self):
        return (
            SpacegroupAnalyzer(self.__bulk).get_conventional_standard_structure()
            if self.use_conventional_cell
            else self.__bulk
        )

    @property
    def generator(self):
        return PymatgenSlabGenerator(
            initial_structure=self.bulk,
            miller_index=self.miller_indices,
            min_slab_size=self.thickness,
            min_vacuum_size=self.vacuum,
            in_unit_planes=True,
            reorient_lattice=True,
        )

    @property
    def __slabs_with_unique_terminations(self):
        return self.generator.get_slabs()

    @property
    def terminations(self):
        return list(set(label_termination(slab) for slab in self.__slabs_with_unique_terminations))

    def get_material(self, termination: str) -> Material:
        for slab in self.__slabs_with_unique_terminations:
            if label_termination(slab) == termination:
                return create_supercell(Material(from_pymatgen(slab)), self.xy_supercell_matrix)
        return None
