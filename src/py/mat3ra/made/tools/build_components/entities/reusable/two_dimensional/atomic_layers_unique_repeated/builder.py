from typing import Type

from ......analyze import BaseMaterialAnalyzer
from ......modify import wrap_to_unit_cell
from ......operations.core.unary import rotate, supercell, translate
from ..... import MaterialWithBuildMetadata
from ....auxiliary.two_dimensional.miller_indices import MillerIndices
from ..crystal_lattice_planes.builder import CrystalLatticePlanesBuilder
from .configuration import AtomicLayersUniqueRepeatedConfiguration


class AtomicLayersUniqueRepeatedBuilder(CrystalLatticePlanesBuilder):
    _ConfigurationType: Type[AtomicLayersUniqueRepeatedConfiguration] = AtomicLayersUniqueRepeatedConfiguration

    def _generate(self, configuration: AtomicLayersUniqueRepeatedConfiguration) -> MaterialWithBuildMetadata:
        crystal_lattice_planes_material = super()._generate(configuration)

        crystal_lattice_planes_material_analyzer = self.get_analyzer(configuration)

        if configuration.termination_top is not None:
            termination = configuration.termination_top
            should_rotate = False
        elif configuration.termination_bottom is not None:
            termination = configuration.termination_bottom
            should_rotate = True
        else:
            raise ValueError(
                "Either termination_top or termination_bottom is required for AtomicLayersUniqueRepeatedBuilder"
            )

        translation_vector = (
            crystal_lattice_planes_material_analyzer.get_translation_vector_for_termination_without_vacuum(termination)
        )
        material_translated = translate(crystal_lattice_planes_material, translation_vector)
        material_translated_wrapped = wrap_to_unit_cell(material_translated)

        if should_rotate:
            # Rotation of basis around [1,0,0] yielded the same material, probably due to the symmetry of unit cell.
            # rotation around X and Z simultaneously gives the mirroring effect. (x,y,z) ⟶ (z,−y,x)
            material_translated_wrapped = rotate(material_translated_wrapped, angle=180, axis=[1, 0, 1])
        material_translated_wrapped_layered = supercell(
            material_translated_wrapped, [[1, 0, 0], [0, 1, 0], [0, 0, configuration.number_of_repetitions]]
        )
        return material_translated_wrapped_layered

    def _update_material_name(
        self, material: MaterialWithBuildMetadata, configuration: AtomicLayersUniqueRepeatedConfiguration
    ) -> MaterialWithBuildMetadata:
        material_analyzer = BaseMaterialAnalyzer(material=material)
        material.formula = material_analyzer.formula
        miller_indices_str = str(MillerIndices(root=configuration.miller_indices))

        if configuration.termination_top is not None:
            termination_str = f"termination {configuration.termination_top}"
        elif configuration.termination_bottom is not None:
            termination_str = f"bottom termination {configuration.termination_bottom}"
        else:
            termination_str = ""

        new_name = f"{material.formula}{miller_indices_str}, {termination_str}"
        material.name = new_name
        return material
