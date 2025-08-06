from typing import Type

from ..... import MaterialWithBuildMetadata

from ......analyze import BaseMaterialAnalyzer
from ....auxiliary.two_dimensional.miller_indices import MillerIndices
from ......modify import wrap_to_unit_cell
from ......operations.core.unary import supercell, translate
from ..crystal_lattice_planes.builder import CrystalLatticePlanesBuilder
from .configuration import AtomicLayersUniqueRepeatedConfiguration


class AtomicLayersUniqueRepeatedBuilder(CrystalLatticePlanesBuilder):
    _ConfigurationType: Type[AtomicLayersUniqueRepeatedConfiguration] = AtomicLayersUniqueRepeatedConfiguration

    def _generate(self, configuration: AtomicLayersUniqueRepeatedConfiguration) -> MaterialWithBuildMetadata:
        crystal_lattice_planes_material = super()._generate(configuration)

        crystal_lattice_planes_material_analyzer = self.get_analyzer(configuration)
        if configuration.termination_top is None:
            raise ValueError("termination_top is required for AtomicLayersUniqueRepeatedBuilder")
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
