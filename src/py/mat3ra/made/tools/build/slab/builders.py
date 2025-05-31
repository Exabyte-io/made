from typing import List, Optional

from mat3ra.esse.models.apse.materials.builders.slab.pymatgen.parameters import (
    PymatgenSlabGeneratorParametersSchema,
)
from pydantic import BaseModel

from mat3ra.made.material import Material
from .configuration import SlabConfiguration, AtomicLayersUniqueRepeated, VacuumConfiguration
from .termination import Termination
from ..supercell import create_supercell
from ..utils import stack_two_components
from ...analyze.other import get_chemical_formula
from ...build import BaseBuilder
from ...build.mixins import ConvertGeneratedItemsPymatgenStructureMixin


class SlabSelectorParameters(BaseModel):
    termination: Termination


class PymatgenSlabGeneratorParameters(PymatgenSlabGeneratorParametersSchema):
    # Parameters described in https://github.com/materialsproject/pymatgen/blob/585bb673c4aa222669c4b0d72ffeec3dbf092630/pymatgen/core/surface.py#L1187
    pass


class SlabBuilderParameters(PymatgenSlabGeneratorParameters):
    make_primitive: bool = False
    use_orthogonal_c: bool = True


class SlabBuilder(ConvertGeneratedItemsPymatgenStructureMixin, BaseBuilder):
    build_parameters: Optional[SlabBuilderParameters]
    _DefaultBuildParameters: SlabBuilderParameters = SlabBuilderParameters()
    _ConfigurationType: type(SlabConfiguration) = SlabConfiguration  # type: ignore
    _GeneratedItemType: Material = Material  # type: ignore
    _SelectorParametersType: type(SlabSelectorParameters) = SlabSelectorParameters  # type: ignore

    def get_material(self, configuration: _ConfigurationType) -> Material:  # type: ignore
        self._configuration = configuration
        return super().get_material(configuration)

    def _generate(self, configuration: _ConfigurationType) -> List[Material]:  # type: ignore
        atomic_layers: AtomicLayersUniqueRepeated = configuration.atomic_layers
        vacuum: VacuumConfiguration = configuration.stack_components[1]
        params = self.build_parameters

        slab_materials = atomic_layers.get_slabs(
            min_slab_size=self._configuration.number_of_layers,
            min_vacuum_size=params.min_vacuum_size,
            in_unit_planes=params.in_unit_planes,
            make_primitive=params.make_primitive,
            symmetrize=params.symmetrize,
            use_orthogonal_c=params.use_orthogonal_c,
        )
        stacked_materials = []
        for slab_material in slab_materials:
            stacked = stack_two_components(slab_material, vacuum, direction=configuration.direction)
            stacked_materials.append(stacked)

        return stacked_materials

    def _post_process(self, items: List[Material], post_process_parameters=None) -> List[Material]:
        materials = items
        materials = [create_supercell(material, self._configuration.xy_supercell_matrix) for material in materials]

        for material in materials:
            if "build" not in material.metadata:
                material.metadata["build"] = {}
            atomic_layers = self._configuration.atomic_layers
            material.metadata["build"]["termination"] = str(atomic_layers.termination_top)
            material.metadata["build"]["configuration"] = self._configuration.to_dict()
            material.metadata["build"]["build_parameters"] = (
                self.build_parameters.model_dump() if self.build_parameters else {}
            )

        return materials

    def _update_material_name(self, material: Material, configuration: SlabConfiguration) -> Material:
        atomic_layers = configuration.atomic_layers

        formula = get_chemical_formula(atomic_layers.crystal)
        miller_indices = "".join(str(i) for i in atomic_layers.miller_indices)
        termination = material.metadata.get("build", {}).get("termination", "")
        material.name = f"{formula}({miller_indices}), termination {termination}, Slab"
        return material
