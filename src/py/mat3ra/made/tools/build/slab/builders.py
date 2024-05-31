from pydantic import BaseModel
from typing import List

from pymatgen.core.surface import SlabGenerator as PymatgenSlabGenerator
from pymatgen.core.interface import label_termination

from mat3ra.made.material import Material
from ...convert import to_pymatgen, PymatgenStructure
from ...build import BaseBuilder
from ...build.mixins import ConvertGeneratedItemsPymatgenStructureMixin
from ..supercell import create_supercell
from .configuration import SlabConfiguration


class SlabSelectorParameters(BaseModel):
    termination: str = ""


class SlabBuildParameters(BaseModel):
    use_orthogonal_z: bool = True
    xy_supercell_matrix: List[List[int]] = [[1, 0], [0, 1]]


class SlabBuilder(BaseBuilder, ConvertGeneratedItemsPymatgenStructureMixin):
    _ConfigurationType: type(SlabConfiguration) = SlabConfiguration  # type: ignore
    _BuildParametersType: type(SlabBuildParameters) = SlabBuildParameters  # type: ignore
    _GeneratedItemType: PymatgenStructure = PymatgenStructure  # type: ignore
    _SelectorParametersType: type(SlabSelectorParameters) = SlabSelectorParameters  # type: ignore

    def _generate(self, configuration: _ConfigurationType) -> List[_GeneratedItemType]:  # type: ignore
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

    def _select(
        self,
        items: List[_GeneratedItemType],
        selector_parameters: _SelectorParametersType = SlabSelectorParameters(),
    ) -> List[_GeneratedItemType]:
        return [slab for slab in items if label_termination(slab) == selector_parameters.termination]

    def _post_process(self, items: List[_GeneratedItemType], post_process_parameters=None) -> List[Material]:
        materials = super()._post_process(items, post_process_parameters)
        return [create_supercell(material, self.build_parameters.xy_supercell_matrix) for material in materials]

    def terminations(self, configuration: _ConfigurationType) -> List[str]:
        return list(set(label_termination(slab) for slab in self._generate_or_get_from_cache(configuration)))
