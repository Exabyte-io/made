from typing import List, Optional, Any, Type

from mat3ra.made.material import Material
from .configurations import (
    CrystalLatticePlanesConfiguration,
    SlabStrainedSupercellConfiguration,
)
from .configurations.base_configurations import AtomicLayersUniqueRepeatedConfiguration
from .configurations.slab_configuration import SlabConfiguration
from .utils import get_orthogonal_c_slab
from .. import BaseBuilderParameters, BaseSingleBuilder, MaterialWithBuildMetadata
from ..stack.builders import StackNComponentsBuilder
from ...analyze.lattice_planes import CrystalLatticePlanesMaterialAnalyzer
from ...analyze.other import get_chemical_formula
from ...modify import wrap_to_unit_cell, translate_to_z_level
from ...operations.core.unary import supercell, translate, strain


class CrystalLatticePlanesBuilder(BaseSingleBuilder):
    _PostProcessParametersType: Any = None
    use_enforce_convention: bool = True

    def _generate(self, configuration: CrystalLatticePlanesConfiguration) -> MaterialWithBuildMetadata:
        crystal_lattice_planes_analyzer = CrystalLatticePlanesMaterialAnalyzer(
            material=configuration.crystal, miller_indices=configuration.miller_indices
        )
        miller_supercell_matrix = crystal_lattice_planes_analyzer.miller_supercell_matrix
        miller_supercell_material = supercell(configuration.crystal, miller_supercell_matrix)
        return miller_supercell_material

    def _enforce_convention(self, material: MaterialWithBuildMetadata) -> MaterialWithBuildMetadata:
        if not self.use_enforce_convention:
            return material
        return translate_to_z_level(material, "bottom")

    def _post_process(
        self, item: Material, post_process_parameters: Optional[_PostProcessParametersType]
    ) -> MaterialWithBuildMetadata:
        item = super()._post_process(item, post_process_parameters)
        return self._enforce_convention(item)


class AtomicLayersUniqueRepeatedBuilder(CrystalLatticePlanesBuilder):
    _ConfigurationType: Type[AtomicLayersUniqueRepeatedConfiguration] = AtomicLayersUniqueRepeatedConfiguration

    def _generate(self, configuration: AtomicLayersUniqueRepeatedConfiguration) -> MaterialWithBuildMetadata:
        crystal_lattice_planes_material = super()._generate(configuration)

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
        return material_translated_wrapped_layered


class SlabBuilderParameters(BaseBuilderParameters):
    use_orthogonal_c: bool = True
    xy_supercell_matrix: List[List[int]] = [[1, 0], [0, 1]]


class SlabBuilder(StackNComponentsBuilder):
    _BuildParametersType = SlabBuilderParameters
    _DefaultBuildParameters: SlabBuilderParameters = SlabBuilderParameters()

    def _configuration_to_material(self, configuration_or_material: Any) -> Material:
        if isinstance(configuration_or_material, AtomicLayersUniqueRepeatedConfiguration):
            builder = AtomicLayersUniqueRepeatedBuilder()
            return builder.get_material(configuration_or_material)
        return super()._configuration_to_material(configuration_or_material)

    def _generate(self, configuration: SlabConfiguration) -> Material:
        stack_as_material = super()._generate(configuration)
        supercell_slab = supercell(stack_as_material, self.build_parameters.xy_supercell_matrix)
        if self.build_parameters.use_orthogonal_c:
            supercell_slab = get_orthogonal_c_slab(supercell_slab)

        return supercell_slab

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
    _ConfigurationType: Type[SlabStrainedSupercellConfiguration] = SlabStrainedSupercellConfiguration

    def _generate(self, configuration: _ConfigurationType) -> Material:
        slab_material = super()._generate(configuration)

        if configuration.xy_supercell_matrix:
            slab_material = supercell(slab_material, configuration.xy_supercell_matrix)
        strained_slab_material = strain(slab_material, configuration.strain_matrix)

        return strained_slab_material
