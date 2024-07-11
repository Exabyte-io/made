from typing import List
from pydantic import BaseModel

from mat3ra.made.material import Material

from ...third_party import PymatgenSlab, PymatgenSlabGenerator, label_pymatgen_slab_termination
from ...analyze import get_chemical_formula
from ...convert import to_pymatgen
from ...build import BaseBuilder
from ...build.mixins import ConvertGeneratedItemsPymatgenStructureMixin
from ..supercell import create_supercell
from .configuration import SlabConfiguration
from .termination import Termination


class SlabSelectorParameters(BaseModel):
    termination: Termination


class SlabBuilder(ConvertGeneratedItemsPymatgenStructureMixin, BaseBuilder):
    _ConfigurationType: type(SlabConfiguration) = SlabConfiguration  # type: ignore
    _GeneratedItemType: PymatgenSlab = PymatgenSlab  # type: ignore
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
            primitive=configuration.make_primitive,
        )
        raw_slabs = generator.get_slabs()
        self.__configuration = configuration

        return [self.__conditionally_convert_slab_to_orthogonal_pymatgen(slab) for slab in raw_slabs]

    def _select(
        self, items: List[_GeneratedItemType], selector_parameters: _SelectorParametersType
    ) -> List[_GeneratedItemType]:
        termination = selector_parameters.termination
        return [slab for slab in items if self.__create_termination_from_slab_pymatgen(slab) == termination]

    def _post_process(self, items: List[_GeneratedItemType], post_process_parameters=None) -> List[Material]:
        materials = super()._post_process(items, post_process_parameters)
        materials = [create_supercell(material, self.__configuration.xy_supercell_matrix) for material in materials]
        for idx, material in enumerate(materials):
            if "build" not in material.metadata:
                material.metadata["build"] = {}
            material.metadata["build"]["termination"] = label_pymatgen_slab_termination(items[idx])

        return materials

    def get_terminations(self, configuration: _ConfigurationType) -> List[Termination]:
        return [
            self.__create_termination_from_slab_pymatgen(slab)
            for slab in self._generate_or_get_from_cache(configuration)
        ]

    def __conditionally_convert_slab_to_orthogonal_pymatgen(self, slab: PymatgenSlab) -> PymatgenSlab:
        return slab.get_orthogonal_c_slab() if self.__configuration.use_orthogonal_z else slab

    def __create_termination_from_slab_pymatgen(self, slab: PymatgenSlab) -> Termination:
        return Termination.from_string(label_pymatgen_slab_termination(slab))

    def _update_material_name(self, material: Material, configuration: SlabConfiguration) -> Material:
        formula = get_chemical_formula(configuration.bulk)
        miller_indices = "".join([str(i) for i in configuration.miller_indices])
        termination = material.metadata.get("build").get("termination", "")
        # for example: "Si8(001), termination Si_P4/mmm_1, Slab"
        new_name = f"{formula}({miller_indices}), termination {termination}, Slab"
        material.name = new_name
        return material
