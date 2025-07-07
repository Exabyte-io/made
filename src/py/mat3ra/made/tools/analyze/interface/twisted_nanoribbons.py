from typing import List, Tuple

import numpy as np
from mat3ra.esse.models.core.abstract.matrix_3x3 import Matrix3x3Schema

from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.interface import InterfaceAnalyzer
from mat3ra.made.tools.analyze.interface.utils.holders import MatchedSubstrateFilmConfigurationHolder
from mat3ra.made.tools.build.slab.configurations import SlabConfiguration
from mat3ra.made.tools.build.slab.configurations import (
    SlabStrainedSupercellConfiguration,
)
from mat3ra.made.tools.modify import rotate, translate_by_vector
from mat3ra.made.utils import get_center_of_coordinates


class TwistedNanoribbonsInterfaceAnalyzer(InterfaceAnalyzer):
    """
    Analyzer for creating twisted interfaces between two nanoribbons.

    Takes two nanoribbons (or any 0D materials) and returns them with:
    - Same lattice vectors
    - One rotated by the specified angle
    - Nanoribbons centered
    - Gap applied using SlabStrainedSupercellWithGapConfiguration

    Args:

        angle (float): Twist angle in degrees.
        vacuum_x (float): Vacuum along x on both sides, in Angstroms.
        vacuum_y (float): Vacuum along y on both sides, in Angstroms.
    """

    angle: float = 0.0
    vacuum_x: float = 5.0
    vacuum_y: float = 5.0

    @property
    def nanoribbon1(self):
        return self.substrate_slab_configuration.atomic_layers.crystal

    @property
    def nanoribbon2(self):
        return self.film_slab_configuration.atomic_layers.crystal

    @property
    def _no_strain_matrix(self) -> Matrix3x3Schema:
        """Return identity matrix for no strain."""
        return Matrix3x3Schema(root=np.eye(3).tolist())

    def _create_primitive_slab_from_nanoribbon(self, nanoribbon: Material) -> SlabConfiguration:
        return SlabConfiguration.from_parameters(
            material_or_dict=nanoribbon,
            miller_indices=(0, 0, 1),
            number_of_layers=1,
            vacuum=0.0,
        )

    def _center_material(self, material: Material) -> Material:
        coordinates = material.basis.coordinates.values
        center_of_mass = get_center_of_coordinates(coordinates)
        lattice_center = np.array([0.5, 0.5, 0.5])
        translation_vector = lattice_center - np.array(center_of_mass)
        return translate_by_vector(material, translation_vector, use_cartesian_coordinates=False)

    def _match_lattice_vectors(self, material1: Material, material2: Material) -> Tuple[Material, Material]:
        lattice1 = material1.lattice.vector_arrays
        lattice2 = material2.lattice.vector_arrays

        n_atoms_1 = len(material1.basis.elements.values)
        n_atoms_2 = len(material2.basis.elements.values)
        if n_atoms_1 >= n_atoms_2:
            target_lattice = lattice1
            material1_matched = material1
            material2_matched = material2.clone()
            material2_matched.set_lattice(material1.lattice)
        else:
            target_lattice = lattice2
            material1_matched = material1.clone()
            material1_matched.set_lattice(material2.lattice)
            material2_matched = material2

        return material1_matched, material2_matched

    def get_strained_configurations(self) -> List[MatchedSubstrateFilmConfigurationHolder]:
        """
        Get strained configurations for the twisted nanoribbon interface.

        Returns:
            List[MatchedSubstrateFilmConfigurationHolder]: List of strained configurations
        """
        centered_nanoribbon1 = self._center_material(self.nanoribbon1)
        centered_nanoribbon2 = self._center_material(self.nanoribbon2)

        matched_nanoribbon1, matched_nanoribbon2 = self._match_lattice_vectors(
            centered_nanoribbon1, centered_nanoribbon2
        )

        rotated_nanoribbon2 = rotate(matched_nanoribbon2, [0, 0, 1], self.angle, wrap=False)

        substrate_slab_config = self._create_primitive_slab_from_nanoribbon(matched_nanoribbon1)
        film_slab_config = self._create_primitive_slab_from_nanoribbon(rotated_nanoribbon2)

        substrate_config = SlabStrainedSupercellConfiguration(
            stack_components=substrate_slab_config.stack_components,
            direction=substrate_slab_config.direction,
            strain_matrix=self._no_strain_matrix,
        )

        film_config = SlabStrainedSupercellConfiguration(
            stack_components=film_slab_config.stack_components,
            direction=film_slab_config.direction,
            strain_matrix=self._no_strain_matrix,
        )

        return [
            MatchedSubstrateFilmConfigurationHolder(
                match_id=0,
                substrate_configuration=substrate_config,
                film_configuration=film_config,
            )
        ]

    @property
    def substrate_nanoribbon_configuration(self):
        return self.get_strained_configurations()[0].substrate_configuration

    @property
    def film_nanoribbon_configuration(self):
        return self.get_strained_configurations()[0].film_configuration
