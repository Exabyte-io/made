from typing import List, Optional

from mat3ra.esse.models.apse.materials.builders.slab.pymatgen.parameters import (
    PymatgenSlabGeneratorParametersSchema,
)
from pydantic import BaseModel

from mat3ra.made.material import Material
from .configuration import SlabConfiguration, AtomicLayersUniqueRepeated, VacuumConfiguration
from .termination import Termination
from ..supercell import create_supercell
from ..utils import stack_two_materials
from ...analyze.other import get_chemical_formula
from ...build import BaseBuilder
from ...build.mixins import ConvertGeneratedItemsPymatgenStructureMixin
from ...convert import to_pymatgen, from_pymatgen
from ...third_party import PymatgenSlab, label_pymatgen_slab_termination


class SlabSelectorParameters(BaseModel):
    termination: Termination


class PymatgenSlabGeneratorParameters(PymatgenSlabGeneratorParametersSchema):
    # Parameters described in https://github.com/materialsproject/pymatgen/blob/585bb673c4aa222669c4b0d72ffeec3dbf092630/pymatgen/core/surface.py#L1187
    pass


class SlabBuilderParameters(PymatgenSlabGeneratorParameters):
    make_primitive: bool = True


class SlabBuilder(ConvertGeneratedItemsPymatgenStructureMixin, BaseBuilder):
    build_parameters: Optional[SlabBuilderParameters] = None
    _ConfigurationType: type(SlabConfiguration) = SlabConfiguration  # type: ignore
    _GeneratedItemType: PymatgenSlab = PymatgenSlab  # type: ignore
    _SelectorParametersType: type(SlabSelectorParameters) = SlabSelectorParameters  # type: ignore

    def get_material(self, configuration: _ConfigurationType) -> Material:  # type: ignore
        self._configuration = configuration
        return super().get_material(configuration)

    @staticmethod
    def _create_vacuum_material(reference: Material, vacuum: VacuumConfiguration) -> Material:
        a_vec, b_vec = reference.lattice.vector_arrays[:2]
        vac_lattice = reference.lattice.from_vectors_array(
            [a_vec, b_vec, [0, 0, vacuum.size]], reference.lattice.units, reference.lattice.type
        )
        return Material.create(
            {
                "name": "Vacuum",
                "lattice": vac_lattice.to_dict(),
                "basis": {"elements": [], "coordinates": []},
                "metadata": {"boundaryConditions": {"type": "pbc", "offset": 0}},
            }
        )

    def _generate(self, configuration: _ConfigurationType) -> List[_GeneratedItemType]:  # type: ignore
        atomic_layers: AtomicLayersUniqueRepeated = configuration.stack_components[0]
        vacuum: VacuumConfiguration = configuration.stack_components[1]
        params = self.build_parameters or SlabBuilderParameters()
        slabs = atomic_layers._generate_pymatgen_slabs(
            min_vacuum_size=params.min_vacuum_size,
            in_unit_planes=True,
            make_primitive=params.make_primitive,
            symmetrize=params.symmetrize,
        )
        if vacuum.is_orthogonal:
            slabs = [slab.get_orthogonal_c_slab() for slab in slabs]
        for idx, slab in enumerate(slabs):
            mat = Material.create(from_pymatgen(slab))
            vacuum_mat = self._create_vacuum_material(mat, vacuum)
            stacked = stack_two_materials(mat, vacuum_mat, direction=configuration.direction)
            slabs[idx] = to_pymatgen(stacked)
        return slabs

    def _post_process(self, items: List[_GeneratedItemType], post_process_parameters=None) -> List[Material]:
        materials = super()._post_process(items, post_process_parameters)
        materials = [create_supercell(material, self._configuration.supercell_xy) for material in materials]

        # Add termination metadata
        for idx, material in enumerate(materials):
            if "build" not in material.metadata:
                material.metadata["build"] = {}
            material.metadata["build"]["termination"] = label_pymatgen_slab_termination(items[idx])
            material.metadata["build"]["configuration"] = self._configuration.to_dict()

        return materials

    def _update_material_name(self, material: Material, configuration: SlabConfiguration) -> Material:
        # Get atomic layers component for material information
        atomic_layers = configuration.stack_components[0]

        formula = get_chemical_formula(atomic_layers.crystal)
        miller_indices = "".join(str(i) for i in atomic_layers.miller_indices)
        termination = material.metadata.get("build", {}).get("termination", "")
        material.name = f"{formula}({miller_indices}), termination {termination}, Slab"
        return material
