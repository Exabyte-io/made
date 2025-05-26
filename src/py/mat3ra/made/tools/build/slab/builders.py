from typing import List, Optional

from mat3ra.made.material import Material
from mat3ra.esse.models.apse.materials.builders.slab.pymatgen.parameters import (
    PymatgenSlabGeneratorParametersSchema,
)
from pydantic import BaseModel

from .configuration import SlabConfiguration, AtomicLayersUniqueRepeated, VacuumConfiguration
from .termination import Termination
from ..supercell import create_supercell
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
    pass


class SlabBuilder(ConvertGeneratedItemsPymatgenStructureMixin, BaseBuilder):
    build_parameters: Optional[SlabBuilderParameters] = None
    _ConfigurationType: type(SlabConfiguration) = SlabConfiguration  # type: ignore
    _GeneratedItemType: PymatgenSlab = PymatgenSlab  # type: ignore
    _SelectorParametersType: type(SlabSelectorParameters) = SlabSelectorParameters  # type: ignore

    def get_material(self, configuration: _ConfigurationType) -> Material:
        self._configuration = configuration
        return super().get_material(configuration)

    def _generate(self, configuration: _ConfigurationType) -> List[_GeneratedItemType]:  # type: ignore
        from ..utils import stack_two_materials

        atomic_layers = configuration.stack_components[0]  # First component is always AtomicLayersUniqueRepeated
        vacuum = configuration.stack_components[1]  # Second component is always VacuumConfiguration

        build_parameters = self.build_parameters or SlabBuilderParameters()
        pymatgen_slabs = atomic_layers._generate_pymatgen_slabs(symmetrize=build_parameters.symmetrize)

        if vacuum.is_orthogonal:
            pymatgen_slabs = [pymatgen_slab.get_orthogonal_c_slab() for pymatgen_slab in pymatgen_slabs]

        # TODO: move to a stacker. Just pass vacuum config
        for i, slab in enumerate(pymatgen_slabs):
            material = Material.create(from_pymatgen(slab))

            # Create vacuum material with slab's in-plane lattice vectors
            slab_vectors = material.lattice.vector_arrays
            vacuum_vectors = [slab_vectors[0], slab_vectors[1], [0, 0, vacuum.size]]
            vacuum_lattice = material.lattice.from_vectors_array(
                vacuum_vectors, material.lattice.units, material.lattice.type
            )
            vacuum_material = Material.create(
                {
                    "name": "Vacuum",
                    "lattice": vacuum_lattice.to_dict(),
                    "basis": {"elements": [], "coordinates": []},
                    "metadata": {"boundaryConditions": {"type": "pbc", "offset": 0}},
                }
            )

            stacked = stack_two_materials(material, vacuum_material, direction=configuration.direction)
            pymatgen_slabs[i] = to_pymatgen(stacked)

        return pymatgen_slabs

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
        miller_indices = "".join([str(i) for i in atomic_layers.miller_indices])
        termination = material.metadata.get("build", {}).get("termination", "")
        # for example: "Si8(001), termination Si_P4/mmm_1, Slab"
        new_name = f"{formula}({miller_indices}), termination {termination}, Slab"
        material.name = new_name
        return material
