from typing import List, Optional, Any, Type

import numpy as np
from mat3ra.esse.models.core.abstract.matrix_3x3 import Matrix3x3Schema

from mat3ra.made.material import Material
from .configurations import (
    CrystalLatticePlanesConfiguration,
    SlabStrainedSupercellConfiguration,
)
from .configurations.base_configurations import AtomicLayersUniqueRepeatedConfiguration
from .configurations.slab_configuration import SlabConfiguration
from .entities import MillerIndices
from .utils import get_orthogonal_c_slab
from .. import BaseBuilderParameters, BaseSingleBuilder, MaterialWithBuildMetadata
from ..stack.builders import StackNComponentsBuilder
from ...analyze import BaseMaterialAnalyzer
from ...analyze.lattice_planes import CrystalLatticePlanesMaterialAnalyzer
from ...modify import wrap_to_unit_cell, translate_to_z_level
from ...operations.core.unary import supercell, translate, strain


class CrystalLatticePlanesBuilder(BaseSingleBuilder):
    _PostProcessParametersType: Any = None
    use_enforce_convention: bool = True

    def get_analyzer(self, configuration: CrystalLatticePlanesConfiguration) -> CrystalLatticePlanesMaterialAnalyzer:
        return CrystalLatticePlanesMaterialAnalyzer(
            material=configuration.crystal, miller_indices=configuration.miller_indices
        )

    def _generate(self, configuration: CrystalLatticePlanesConfiguration) -> MaterialWithBuildMetadata:
        crystal_lattice_planes_analyzer = self.get_analyzer(configuration)
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

        crystal_lattice_planes_material_analyzer = self.get_analyzer(configuration)
        translation_vector = (
            crystal_lattice_planes_material_analyzer.get_translation_vector_for_termination_without_vacuum(
                configuration.termination_top
            )
        )
        material_translated = translate(crystal_lattice_planes_material, translation_vector)
        material_translated_wrapped = wrap_to_unit_cell(material_translated)
        material_translated_wrapped_layered = supercell(
            material_translated_wrapped, [[1, 0, 0], [0, 1, 0], [0, 0, configuration.number_of_repetitions]]
        )
        return material_translated_wrapped_layered

    def _update_material_name(
        self, material: MaterialWithBuildMetadata, configuration: AtomicLayersUniqueRepeatedConfiguration
    ) -> MaterialWithBuildMetadata:
        material_analyzer = BaseMaterialAnalyzer(material=material)
        material.formula = material_analyzer.formula
        termination = configuration.termination_top
        miller_indices_str = str(MillerIndices(root=configuration.miller_indices))
        # for example: "Si(001), termination Si_P4/mmm_1"
        new_name = f"{material.formula}{miller_indices_str}, termination {termination}"
        material.name = new_name
        return material


class SlabBuilderParameters(BaseBuilderParameters):
    use_orthogonal_c: bool = True
    xy_supercell_matrix: List[List[int]] = [[1, 0], [0, 1]]


class SlabBuilder(StackNComponentsBuilder):
    _BuildParametersType = SlabBuilderParameters
    _DefaultBuildParameters: SlabBuilderParameters = SlabBuilderParameters()

    @property
    def stack_component_types_conversion_map(self):
        return {
            **super().stack_component_types_conversion_map,
            AtomicLayersUniqueRepeatedConfiguration: AtomicLayersUniqueRepeatedBuilder,
        }

    def _generate(self, configuration: SlabConfiguration) -> Material:
        stack_as_material = super()._generate(configuration)
        supercell_slab = supercell(stack_as_material, self.build_parameters.xy_supercell_matrix)
        if self.build_parameters.use_orthogonal_c:
            supercell_slab = get_orthogonal_c_slab(supercell_slab)
        return supercell_slab

    def _update_material_name(self, material: MaterialWithBuildMetadata, configuration: SlabConfiguration) -> Material:
        # for example: "Si(001), termination Si_P4/mmm_1, Slab"
        material = AtomicLayersUniqueRepeatedBuilder()._update_material_name(material, configuration.atomic_layers)
        new_name = f"{material.name}, Slab"
        material.name = new_name
        return material


class SlabStrainedSupercellBuilder(SlabBuilder):
    _ConfigurationType: Type[SlabStrainedSupercellConfiguration] = SlabStrainedSupercellConfiguration

    def _generate(self, configuration: _ConfigurationType) -> Material:
        slab_material = super()._generate(configuration)
        # TODO: Fix pymatgen's strain matching (issue for HEX a and b switching)
        diagonal_strain_matrix = [
            configuration.strain_matrix.root[i].root[i] for i in range(len(configuration.strain_matrix.root))
        ]

        if configuration.xy_supercell_matrix:
            slab_material = supercell(slab_material, configuration.xy_supercell_matrix)

        # In case pymatgen returned strain matrix that has supercell matrix embedded, we extract it.
        integer_parts = [int(val) if val >= 1 else 1 for val in diagonal_strain_matrix]
        fractional_parts = [val / int(val) if val >= 1 else val for val in diagonal_strain_matrix]

        # Create supercell if any integer parts are > 1
        if any(val > 1 for val in integer_parts):
            supercell_matrix = np.diag(integer_parts).tolist()
            slab_material = supercell(slab_material, supercell_matrix)

        # Apply remaining fractional strain
        diagonal_strain_3x3 = np.diag(fractional_parts)
        diagonal_strain_schema = Matrix3x3Schema(root=diagonal_strain_3x3.tolist())
        strained_slab_material = strain(slab_material, diagonal_strain_schema)

        return strained_slab_material
