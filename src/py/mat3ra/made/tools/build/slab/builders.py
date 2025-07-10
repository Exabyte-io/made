from typing import List, Optional, Any, Type

import numpy as np
from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum

from mat3ra.made.material import Material
from .configurations import (
    CrystalLatticePlanesConfiguration,
    SlabStrainedSupercellConfiguration,
    SlabStrainedSupercellWithGapConfiguration,
)
from .configurations.base_configurations import AtomicLayersUniqueRepeatedConfiguration
from .configurations.slab_configuration import SlabConfiguration
from .utils import get_orthogonal_c_slab
from .. import BaseBuilderParameters, BaseSingleBuilder
from ..stack.builders import Stack2ComponentsBuilder
from ...analyze.lattice_planes import CrystalLatticePlanesMaterialAnalyzer
from ...analyze.other import get_chemical_formula, get_atomic_coordinates_extremum
from ...modify import wrap_to_unit_cell, translate_to_z_level
from ...operations.core.unary import supercell, translate, strain, edit_cell
from ...utils import AXIS_TO_INDEX_MAP


class CrystalLatticePlanesBuilder(BaseSingleBuilder):
    _PostProcessParametersType: Any = None
    use_enforce_convention: bool = True

    def _generate(self, configuration: CrystalLatticePlanesConfiguration) -> Material:
        crystal_lattice_planes_analyzer = CrystalLatticePlanesMaterialAnalyzer(
            material=configuration.crystal, miller_indices=configuration.miller_indices
        )
        miller_supercell_matrix = crystal_lattice_planes_analyzer.miller_supercell_matrix
        miller_supercell_material = supercell(configuration.crystal, miller_supercell_matrix)
        return miller_supercell_material

    def _enforce_convention(self, material: Material) -> Material:
        if not self.use_enforce_convention:
            return material
        return translate_to_z_level(material, "bottom")

    def _post_process(self, item: Material, post_process_parameters: Optional[_PostProcessParametersType]) -> Material:
        item = super()._post_process(item, post_process_parameters)
        return self._enforce_convention(item)


class AtomicLayersUniqueRepeatedBuilder(CrystalLatticePlanesBuilder):
    _ConfigurationType: Type[AtomicLayersUniqueRepeatedConfiguration] = AtomicLayersUniqueRepeatedConfiguration

    def _generate(self, configuration: AtomicLayersUniqueRepeatedConfiguration) -> Material:
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


class SlabBuilder(Stack2ComponentsBuilder):
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


class SlabWithGapBuilder(SlabStrainedSupercellBuilder):
    _ConfigurationType: Type[SlabStrainedSupercellWithGapConfiguration] = SlabStrainedSupercellWithGapConfiguration

    def _generate(self, configuration: _ConfigurationType) -> Material:
        strained_slab_material = super()._generate(configuration)

        if configuration.gap is not None:
            strained_slab_material = self._adjust_lattice_for_gap(
                strained_slab_material, configuration.gap, configuration.gap_direction
            )

        return strained_slab_material

    def _adjust_lattice_for_gap(self, material: Material, gap: float, direction: AxisEnum) -> Material:
        """
        Adjust the cell along the stacking direction to make the distance from the cell end to its closest atom
        to be equal to the gap, in Angstroms.
        """
        direction_str = direction.value
        axis_index = AXIS_TO_INDEX_MAP[direction_str]

        max_fractional = get_atomic_coordinates_extremum(material, "max", direction_str, False)
        current_vectors = material.lattice.vector_arrays
        current_vector = np.array(current_vectors[axis_index])
        current_length = np.linalg.norm(current_vector)

        new_length = (max_fractional * current_length) + gap

        if current_length > 0:
            new_vector = (current_vector / current_length) * new_length
        else:
            new_vector = current_vector
        new_lattice_vectors = list(current_vectors)
        new_lattice_vectors[axis_index] = new_vector.tolist()
        new_material = edit_cell(material, new_lattice_vectors)
        return new_material
