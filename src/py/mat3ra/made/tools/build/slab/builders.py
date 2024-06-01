from pymatgen.core.surface import SlabGenerator as PymatgenSlabGenerator
from pymatgen.core.interface import label_termination
from typing import List
from pydantic import BaseModel


from mat3ra.made.material import Material

from ...analyze import get_chemical_formula
from ...convert import to_pymatgen, PymatgenStructure, PymatgenSlab
from ...build import BaseBuilder
from ...build.mixins import ConvertGeneratedItemsPymatgenStructureMixin
from ..supercell import create_supercell
from .configuration import SlabConfiguration


class SlabSelectorParameters(BaseModel):
    termination: str = ""


class SlabBuilder(ConvertGeneratedItemsPymatgenStructureMixin, BaseBuilder):
    _ConfigurationType: type(SlabConfiguration) = SlabConfiguration  # type: ignore
    _GeneratedItemType: PymatgenStructure = PymatgenStructure  # type: ignore
    _SelectorParametersType: type(SlabSelectorParameters) = SlabSelectorParameters  # type: ignore
    __configuration: SlabConfiguration

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
        self.__configuration = configuration

        return [self.__conditionally_convert_slab_to_orthogonal_pymatgen(slab) for slab in raw_slabs]

    def _select(
        self,
        items: List[_GeneratedItemType],
        selector_parameters: _SelectorParametersType = SlabSelectorParameters(),
    ) -> List[_GeneratedItemType]:
        return [slab for slab in items if label_termination(slab) == selector_parameters.termination]

    def _post_process(self, items: List[_GeneratedItemType], post_process_parameters=None) -> List[Material]:
        materials = super()._post_process(items, post_process_parameters)
        materials = [create_supercell(material, self.__configuration.xy_supercell_matrix) for material in materials]
        for idx, material in enumerate(materials):
            material.metadata["termination"] = label_termination(items[idx])

        return materials

    def get_terminations(self, configuration: _ConfigurationType) -> List[str]:
        return list(set(label_termination(slab) for slab in self._generate_or_get_from_cache(configuration)))

    def __conditionally_convert_slab_to_orthogonal_pymatgen(self, slab: PymatgenSlab) -> PymatgenStructure:
        return slab.get_orthogonal_c_slab() if self.__configuration.use_orthogonal_z else slab

    def _update_material_name(self, material: Material, configuration: SlabConfiguration) -> Material:
        formula = get_chemical_formula(configuration.bulk)
        miller_indices = "".join([str(i) for i in configuration.miller_indices])
        termination = material.metadata["termination"]
        # for example: "Si8(001), termination Si_P4/mmm_1, Slab"
        new_name = f"{formula}({miller_indices}), termination {termination}, Slab"
        material.name = new_name
        return material
