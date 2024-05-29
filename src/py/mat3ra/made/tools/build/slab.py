from pymatgen.core.structure import Structure
from pymatgen.core.surface import SlabGenerator as PymatgenSlabGenerator
from pymatgen.core.interface import label_termination
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from typing import List, Tuple, Optional

from . import BaseBuilder
from .supercell import create_supercell
from ..convert import from_pymatgen, to_pymatgen
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
        self.__bulk: Structure = (
            SpacegroupAnalyzer(to_pymatgen(bulk)).get_conventional_standard_structure()
            if self.use_conventional_cell
            else to_pymatgen(bulk)
        )
        self.__miller_indices = miller_indices
        self.__builder = SlabBuilder()
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
    def terminations(self):
        return self.__builder.terminations(self)

    def get_material(self, termination: Optional[str]):
        return self.__builder.get_material(self, termination)


class SlabBuilder(BaseBuilder):
    def __init__(self):
        pass

    def get_material(
        self, configuration: SlabConfiguration = SlabConfiguration(), termination: Optional[str] = "", **kwargs
    ) -> Material:
        for slab in SlabBuilder.__slabs_with_unique_terminations(self, configuration):
            if label_termination(slab) == termination:
                return create_supercell(
                    Material(from_pymatgen(slab)),
                    configuration.xy_supercell_matrix,
                )
        raise ValueError(f"Termination {termination} not found in slabs.")

    def __generator(self, slab_configuration: SlabConfiguration) -> PymatgenSlabGenerator:
        return PymatgenSlabGenerator(
            initial_structure=to_pymatgen(slab_configuration.bulk),
            miller_index=slab_configuration.miller_indices,
            min_slab_size=slab_configuration.thickness,
            min_vacuum_size=slab_configuration.vacuum,
            in_unit_planes=True,
            reorient_lattice=True,
        )

    def __slabs_with_unique_terminations(self, configuration: SlabConfiguration):
        return [
            slab.get_orthogonal_c_slab() if configuration.use_orthogonal_z else slab
            for slab in self.__generator(configuration).get_slabs()
        ]

    def terminations(self, slab_configuration: SlabConfiguration) -> List[str]:
        return list(set(label_termination(slab) for slab in self.__slabs_with_unique_terminations(slab_configuration)))
