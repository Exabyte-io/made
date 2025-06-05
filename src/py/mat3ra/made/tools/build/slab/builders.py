from mat3ra.code.vector import Vector3D

from mat3ra.made.material import Material
from .configuration import (
    SlabConfiguration,
    AtomicLayersUniqueRepeatedConfiguration,
    CrystalLatticePlanesConfiguration,
    ConventionalCellConfiguration,
)
from ..stack.builders import StackBuilder2Components
from ...build import BaseBuilder
from ...convert import to_pymatgen, from_pymatgen
from ...operations.core.unary import supercell, translate
from ...third_party import PymatgenSpacegroupAnalyzer


class ConventionalCellBuilder(BaseBuilder):
    def get_material(self, configuration: ConventionalCellConfiguration) -> Material:
        bulk_pmg = to_pymatgen(configuration.crystal)
        conventional_cell = PymatgenSpacegroupAnalyzer(bulk_pmg).get_conventional_standard_structure()
        return Material.create(from_pymatgen(conventional_cell))


class CrystalLatticePlanesBuilder(BaseBuilder):
    @staticmethod
    def calculate_orthogonal_lattice(crystal, miller_supercell):
        return supercell(crystal, miller_supercell)

    def get_material(self, configuration: CrystalLatticePlanesConfiguration) -> Material:
        orthogonal_c_cell = self.calculate_orthogonal_lattice(configuration.crystal, configuration.miller_supercell)
        return supercell(orthogonal_c_cell, configuration.rotational_matrix)


class AtomicLayersUniqueRepeatedBuilder(BaseBuilder):
    _ConfigurationType = AtomicLayersUniqueRepeatedConfiguration

    def get_material(self, configuration: _ConfigurationType) -> Material:
        material = configuration.orthogonal_c_cell
        translation_vector: Vector3D = configuration.get_translation_vector(configuration.termination_top)
        material_translated = translate(material, translation_vector)
        material_translated_with_repetitions = supercell(
            material_translated, [[1, 0, 0], [0, 1, 0], [0, 0, configuration.number_of_repetitions]]
        )
        return material_translated_with_repetitions


class SlabBuilderParameters(BaseBuilder):
    make_primitive: bool = False
    use_orthogonal_c: bool = True


class SlabBuilder(StackBuilder2Components):

    def generate(self, configuration: SlabConfiguration) -> Material:
        stacked_materials = super().generate(configuration)
        supercell_slab = supercell(stacked_materials, configuration.xy_supercell_matrix)
        return supercell_slab

    def get_material(self, configuration: SlabConfiguration) -> Material:
        return self.generate(configuration)
