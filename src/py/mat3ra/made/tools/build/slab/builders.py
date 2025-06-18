from typing import List, Optional, Any, Type

from mat3ra.made.material import Material
from .configuration import (
    SlabConfiguration,
    AtomicLayersUniqueRepeatedConfiguration,
    CrystalLatticePlanesConfiguration,
)
from .utils import get_orthogonal_c_slab
from .. import BaseBuilderParameters
from ..slab.configuration import SlabStrainedSupercellConfiguration
from ..stack.builders import StackBuilder2Components
from ...analyze.lattice_planes import CrystalLatticePlanesMaterialAnalyzer
from ...analyze.other import get_chemical_formula
from ...build import BaseBuilder
from ...modify import wrap_to_unit_cell, translate_to_z_level
from ...operations.core.unary import supercell, translate, strain


class CrystalLatticePlanesBuilder(BaseBuilder):
    _GeneratedItemType: Material = Material
    _PostProcessParametersType: Any = None
    use_enforce_convention: bool = True

    def _generate(self, configuration: CrystalLatticePlanesConfiguration) -> List[Material]:
        crystal_lattice_planes_analyzer = CrystalLatticePlanesMaterialAnalyzer(
            material=configuration.crystal, miller_indices=configuration.miller_indices
        )
        miller_supercell_matrix = crystal_lattice_planes_analyzer.miller_supercell_matrix
        miller_supercell_material = supercell(configuration.crystal, miller_supercell_matrix)
        return [miller_supercell_material]

    def _enforce_convention(self, material: Material) -> Material:
        if not self.use_enforce_convention:
            return material
        return translate_to_z_level(material, "bottom")

    def _post_process(
        self, items: List[_GeneratedItemType], post_process_parameters: Optional[_PostProcessParametersType]
    ) -> List[Material]:
        items = super()._post_process(items, post_process_parameters)
        return [self._enforce_convention(i) for i in items]


class AtomicLayersUniqueRepeatedBuilder(CrystalLatticePlanesBuilder):
    _ConfigurationType: Type[AtomicLayersUniqueRepeatedConfiguration] = AtomicLayersUniqueRepeatedConfiguration

    def _generate(self, configuration: AtomicLayersUniqueRepeatedConfiguration) -> List[Material]:
        crystal_lattice_planes_material = super()._generate(configuration)[0]

        crystal_lattice_planes_analyzer = CrystalLatticePlanesMaterialAnalyzer(
            material=configuration.crystal, miller_indices=configuration.miller_indices
        )
        translation_vector = crystal_lattice_planes_analyzer.get_translation_vector_for_termination_without_vacuum(
            configuration.termination_top
        )
        material_translated = translate(crystal_lattice_planes_material, translation_vector)
        material_translated_wrapped = wrap_to_unit_cell(material_translated)
        material_translated_wrapped_layered = supercell(
            material_translated_wrapped, [[1, 0, 0], [0, 1, 0], [0, 0, configuration.number_of_repetitions]]
        )
        return [material_translated_wrapped_layered]


class SlabBuilderParameters(BaseBuilderParameters):
    use_orthogonal_c: bool = True
    xy_supercell_matrix: List[List[int]] = [[1, 0], [0, 1]]


class SlabBuilder(StackBuilder2Components):
    _BuildParametersType = SlabBuilderParameters
    _DefaultBuildParameters: SlabBuilderParameters = SlabBuilderParameters()

    def _generate(self, configuration: SlabConfiguration) -> List[Material]:
        stack_as_material_list = super()._generate(configuration)
        stack_as_material = stack_as_material_list[0]
        supercell_slab = supercell(stack_as_material, self.build_parameters.xy_supercell_matrix)
        if self.build_parameters.use_orthogonal_c:
            supercell_slab = get_orthogonal_c_slab(supercell_slab)

        return [supercell_slab]

    def _update_material_name(self, material: Material, configuration: SlabConfiguration) -> Material:
        atomic_layers = configuration.atomic_layers

        formula = get_chemical_formula(configuration.atomic_layers.crystal)
        miller_indices_str = "".join([str(i) for i in atomic_layers.miller_indices])
        termination = atomic_layers.termination_top

        # for example: "Si8(001), termination Si_P4/mmm_1, Slab"
        new_name = f"{formula}({miller_indices_str}), termination {termination}, Slab"
        material.name = new_name
        return material


class SlabStrainedSupercellBuilder(SlabBuilder):
    def _generate(self, configuration: SlabStrainedSupercellConfiguration) -> List[Material]:
        materials = super()._generate(configuration)

        material = materials[0]
        if configuration.xy_supercell_matrix:
            material = supercell(material, configuration.xy_supercell_matrix)
        if configuration.strain_matrix:
            material = strain(material, configuration.strain_matrix)

        return [material]
