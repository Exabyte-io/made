from pymatgen.core.structure import Structure as PymatgenStructure
from pymatgen.core.surface import SlabGenerator as PymatgenSlabGenerator
from pymatgen.core.interface import label_termination
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer as PymatgenSpacegroupAnalyzer
from typing import List, Tuple
from pydantic import BaseModel

from . import BaseBuilder
from .supercell import create_supercell
from ..convert import from_pymatgen, to_pymatgen
from ...material import Material


class BaseSlabConfiguration(BaseModel):
    @property
    def bulk(self) -> Material:
        raise NotImplementedError

    @property
    def miller_indices(self) -> Tuple[int, int, int]:
        raise NotImplementedError


class SlabBuildParameters(BaseModel):
    use_orthogonal_z: bool = True
    xy_supercell_matrix: List[List[int]] = [[1, 0], [0, 1]]


class SlabSelectorParameters(BaseModel):
    termination: str = ""


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


class SlabBuilder(BaseBuilder[SlabConfiguration, SlabBuildParameters, PymatgenStructure, SlabSelectorParameters, None]):

    def __generate(self, configuration: SlabConfiguration) -> List[PymatgenStructure]:
        generator = PymatgenSlabGenerator(
            initial_structure=to_pymatgen(configuration.bulk),
            miller_index=configuration.miller_indices,
            min_slab_size=configuration.thickness,
            min_vacuum_size=configuration.vacuum,
            in_unit_planes=True,
            reorient_lattice=True,
        )
        raw_slabs = generator.get_slabs()

        slabs = [slab.get_orthogonal_c_slab() if self.build_parameters.use_orthogonal_z else slab for slab in raw_slabs]

        return slabs

    def __select(
        self,
        items: List[PymatgenStructure],
        selector_parameters: SlabSelectorParameters = SlabSelectorParameters(),
    ) -> List[PymatgenStructure]:
        return [slab for slab in items if label_termination(slab) == selector_parameters.termination]

    def __post_process(self, items: List[PymatgenStructure], post_process_parameters=None) -> List[Material]:
        materials = [Material(from_pymatgen(slab)) for slab in items]
        return [create_supercell(material, self.build_parameters.xy_supercell_matrix) for material in materials]

    def terminations(self, configuration: SlabConfiguration) -> List[str]:
        # TODO: use __generate_or_get_from_cache()
        return list(set(label_termination(slab) for slab in self.__generate(configuration)))
