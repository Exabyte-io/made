from typing import List, Optional

from mat3ra.made.material import Material
from mat3ra.esse.models.apse.materials.builders.slab.pymatgen.parameters import (
    PymatgenSlabGeneratorParametersSchema,
)
from pydantic import BaseModel

from .configuration import SlabConfiguration, CrystalLatticePlanes, AtomicLayersUniqueRepeated, VacuumConfiguration
from .termination import Termination
from ..supercell import create_supercell
from ..utils import stack_two_components
from ...analyze.other import get_chemical_formula
from ...build import BaseBuilder
from ...build.mixins import ConvertGeneratedItemsPymatgenStructureMixin
from ...convert import to_pymatgen
from ...modify import add_vacuum
from ...third_party import PymatgenSlab, PymatgenSlabGenerator, label_pymatgen_slab_termination


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
        """
        Get a material from the configuration.

        Args:
            configuration (SlabConfiguration): The configuration to use.

        Returns:
            Material: The generated material.
        """
        self._configuration = configuration
        return super().get_material(configuration)

    def _generate(self, configuration: _ConfigurationType) -> List[_GeneratedItemType]:  # type: ignore
        """
        Generate a slab by stacking components with proper vacuum and terminations.

        Args:
            configuration (SlabConfiguration): The configuration containing stack components and parameters.

        Returns:
            List[PymatgenSlab]: List of generated slabs.
        """
        # Stack all components in order
        components = configuration.stack_components
        stacked = components[0]
        for comp in components[1:]:
            stacked = stack_two_components(stacked, comp, direction=configuration.direction)
        # If the result is not a PymatgenSlab, convert it if needed
        pymatgen_structure = to_pymatgen(stacked)
        return [pymatgen_structure]

    def _post_process(self, items: List[_GeneratedItemType], post_process_parameters=None) -> List[Material]:
        materials = super()._post_process(items, post_process_parameters)
        materials = [create_supercell(material, self._configuration.supercell_xy) for material in materials]
        build_parameters = self.build_parameters or SlabBuilderParameters()

        # Get components from stack components
        atomic_layers = next(component for component in self._configuration.stack_components 
                           if isinstance(component, AtomicLayersUniqueRepeated))
        vacuum = next(component for component in self._configuration.stack_components 
                     if isinstance(component, VacuumConfiguration))

        # Adding total vacuum to be exactly as specified in configuration
        added_vacuum = (
            build_parameters.min_vacuum_size * atomic_layers.crystal.lattice.c
            if build_parameters.in_unit_planes
            else build_parameters.min_vacuum_size
        )
        vacuum_to_add = vacuum.size - added_vacuum

        materials_with_vacuum = [add_vacuum(material, vacuum_to_add) for material in materials]
        for idx, material in enumerate(materials_with_vacuum):
            if "build" not in material.metadata:
                material.metadata["build"] = {}
            material.metadata["build"]["termination"] = label_pymatgen_slab_termination(items[idx])

        return materials_with_vacuum

    def _update_material_name(self, material: Material, configuration: SlabConfiguration) -> Material:
        # Get atomic layers component for material information
        atomic_layers = next(component for component in configuration.stack_components 
                           if isinstance(component, AtomicLayersUniqueRepeated))
        
        formula = get_chemical_formula(atomic_layers.crystal)
        miller_indices = "".join([str(i) for i in atomic_layers.miller_indices])
        termination = material.metadata.get("build", {}).get("termination", "")
        # for example: "Si8(001), termination Si_P4/mmm_1, Slab"
        new_name = f"{formula}({miller_indices}), termination {termination}, Slab"
        material.name = new_name
        return material
