from typing import List

import numpy as np
from mat3ra.code.vector import Vector3D

from mat3ra.made.material import Material
from .configuration import (
    SlabConfiguration,
    AtomicLayersUniqueRepeatedConfiguration,
    CrystalLatticePlanesConfiguration,
    ConventionalCellConfiguration,
    TerminationAnalyzer,
)
from .. import BaseBuilderParameters
from ..stack.builders import StackBuilder2Components
from ...analyze.other import get_chemical_formula
from ...build import BaseBuilder
from ...convert import to_pymatgen, from_pymatgen
from ...modify import wrap_to_unit_cell
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

    def get_surface_supercell(self, configuration: _ConfigurationType) -> Material:
        return supercell(configuration.crystal, configuration.miller_supercell)

    def get_translation_vector(self, configuration) -> Vector3D:
        termination_analyzer = TerminationAnalyzer(configuration.crystal, configuration.miller_indices)
        vector_crystal = termination_analyzer.find_translation_vector(configuration.termination_top)
        vector_cartesian = self.get_surface_supercell(configuration).basis.cell.convert_point_to_cartesian(
            vector_crystal
        )
        return vector_cartesian

    def get_material(self, configuration: _ConfigurationType) -> Material:
        material = self.get_surface_supercell(configuration)
        translation_vector: Vector3D = self.get_translation_vector(configuration)
        material_translated = translate(material, translation_vector)
        material_translated_wrapped = wrap_to_unit_cell(material_translated)
        material_translated_with_repetitions = supercell(
            material_translated_wrapped, [[1, 0, 0], [0, 1, 0], [0, 0, configuration.number_of_repetitions]]
        )
        return material_translated_with_repetitions


class SlabBuilderParameters(BaseBuilderParameters):
    use_orthogonal_c: bool = True
    xy_supercell_matrix: List[List[int]] = [[1, 0], [0, 1]]


class SlabBuilder(StackBuilder2Components):
    _BuildParametersType = SlabBuilderParameters
    DefaultBuilderParameters: SlabBuilderParameters = SlabBuilderParameters()

    def _generate(self, configuration: SlabConfiguration) -> List[Material]:
        stack_as_material_list = super()._generate(configuration)
        stack_as_material = stack_as_material_list[0]

        supercell_slab = supercell(stack_as_material, self.build_parameters.xy_supercell_matrix)
        if self.build_parameters.use_orthogonal_c:
            supercell_slab = self._make_orthogonal_c(supercell_slab)
        return [supercell_slab]

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

    def _update_material_name(self, material: Material, configuration: SlabConfiguration) -> Material:
        atomic_layers = configuration.atomic_layers

        formula = get_chemical_formula(configuration.atomic_layers.crystal)
        miller_indices_str = "".join([str(i) for i in atomic_layers.miller_indices])
        termination = atomic_layers.termination_top

        new_name = f"{formula}({miller_indices_str}), termination {termination}, Slab"
        material.name = new_name
        return material
