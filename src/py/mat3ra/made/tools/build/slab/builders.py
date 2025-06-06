import numpy as np
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
from ...operations.core.unary import supercell, translate, edit_cell
from ...third_party import PymatgenSpacegroupAnalyzer


class ConventionalCellBuilder(BaseBuilder):
    def get_material(self, configuration: ConventionalCellConfiguration) -> Material:
        bulk_pmg = to_pymatgen(configuration.crystal)
        conventional_cell = PymatgenSpacegroupAnalyzer(bulk_pmg).get_conventional_standard_structure()
        return Material.create(from_pymatgen(conventional_cell))


class CrystalLatticePlanesBuilder(BaseBuilder):
    @staticmethod
    def calculate_plane_lattice(crystal, miller_supercell):
        return supercell(crystal, miller_supercell)

    def get_material(self, configuration: CrystalLatticePlanesConfiguration) -> Material:
        plane_cell = self.calculate_plane_lattice(configuration.crystal, configuration.miller_supercell)
        # return supercell(plane_cell, configuration.rotational_matrix)
        return plane_cell


class AtomicLayersUniqueRepeatedBuilder(BaseBuilder):
    _ConfigurationType = AtomicLayersUniqueRepeatedConfiguration

    def get_material(self, configuration: _ConfigurationType) -> Material:
        material = configuration.surface_supercell
        translation_vector: Vector3D = configuration.get_translation_vector(configuration.termination_top)
        material_translated = translate(material, translation_vector)
        material_translated_with_repetitions = supercell(
            material_translated, [[1, 0, 0], [0, 1, 0], [0, 0, configuration.number_of_repetitions]]
        )
        return material_translated_with_repetitions


class SlabBuilder(StackBuilder2Components):

    def generate(self, configuration: SlabConfiguration) -> Material:
        stacked_materials = super().generate(configuration)
        supercell_slab = supercell(stacked_materials, configuration.xy_supercell_matrix)
        if configuration.use_orthogonal_c:
            supercell_slab = self._make_orthogonal_c(supercell_slab)
        return supercell_slab

    def get_material(self, configuration: SlabConfiguration) -> Material:
        return self.generate(configuration)

    def _make_orthogonal_c(self, material: Material) -> Material:
        """
        Make the c-vector orthogonal to the ab plane after slab construction is complete.
        This should be applied after vacuum has been added to the slab.
        """
        current_vectors = material.lattice.vector_arrays
        a = np.array(current_vectors[0])
        b = np.array(current_vectors[1])
        c_old = np.array(current_vectors[2])

        normal = np.cross(a, b)
        norm_norm = np.linalg.norm(normal)
        if norm_norm < 1e-8:
            raise ValueError("Vectors a and b are collinear or too small to define a plane.")
        n_hat = normal / norm_norm

        height = float(np.dot(c_old, n_hat))

        new_orthogonal_vector_c = (n_hat * height).tolist()

        new_vectors = [
            current_vectors[0],
            current_vectors[1],
            new_orthogonal_vector_c,
        ]

        return edit_cell(material, new_vectors)
