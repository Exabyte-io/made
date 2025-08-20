from functools import cached_property
from typing import List, Optional

import numpy as np
from mat3ra.code.entity import InMemoryEntityPydantic
from mat3ra.esse.models.core.abstract.matrix_3x3 import Matrix3x3Schema
from mat3ra.esse.models.materials_category_components.entities.auxiliary.two_dimensional.supercell_matrix_2d import (
    SupercellMatrix2DSchema,
)
from mat3ra.made.lattice import Lattice
from mat3ra.made.material import Material
from mat3ra.utils.matrix import convert_2x2_to_3x3

from ....utils import create_2d_supercell_matrices
from ..utils import calculate_von_mises_strain
from .simple import InterfaceAnalyzer
from .utils.holders import MatchedSubstrateFilmConfigurationHolder
from .utils.vector import align_first_vector_to_x_2d_right_handed, are_vectors_colinear


class CSLMatchHolder(InMemoryEntityPydantic):
    match_id: int
    substrate_supercell_matrix: SupercellMatrix2DSchema
    film_supercell_matrix: SupercellMatrix2DSchema
    film_strain_matrix: Matrix3x3Schema
    match_area: float
    total_strain_percentage: float
    film_rotation_angle: float


class CSLInterfaceAnalyzer(InterfaceAnalyzer):
    """
    Interface analyzer using Coincidence Site Lattice (CSL) method.

    This analyzer combines aspects of ZSL and commensurate methods by:
    1. Creating slabs for both substrate and film
    2. Ensuring both have the same lattice type
    3. Aligning lattice vectors (a_sub || a_film, b_sub || b_film)
    4. Applying rotational supercells to film at various angles
    5. Finding substrate diagonal matrices that match film within tolerance
    6. Calculating strain and creating matched configurations

    Attributes:
        max_area (float): Maximum area for supercell matching
        length_tolerance (float): Tolerance for matching lattice vector lengths
        angle_step (float): Step size for rotation angles in degrees
        max_rotation_angle (float): Maximum rotation angle to test in degrees
        max_supercell_size (int): Maximum supercell matrix element size
        strain_tolerance (float): Maximum acceptable strain percentage
    """

    max_area: float = 100.0
    length_tolerance: float = 0.05
    angle_step: float = 5.0
    max_rotation_angle: float = 180.0
    max_supercell_size: int = 10
    strain_tolerance: float = 10.0

    def _validate_same_lattice_type(self) -> bool:
        """Check if substrate and film have the same lattice type."""
        substrate_type = self.substrate_material.lattice.type
        film_type = self.film_material.lattice.type
        return substrate_type == film_type

    def _get_lattice_vectors_2d(self, material: Material) -> np.ndarray:
        """Get 2D lattice vectors (a, b) from material."""
        return np.array(material.lattice.vector_arrays[:2])[:, :2]

    def _create_rotation_matrix(self, angle_degrees: float) -> np.ndarray:
        """Create 2D rotation matrix for given angle in degrees."""
        angle_rad = np.radians(angle_degrees)
        cos_a = np.cos(angle_rad)
        sin_a = np.sin(angle_rad)
        return np.array([[cos_a, -sin_a], [sin_a, cos_a]])

    def _apply_supercell_and_rotation(
        self, vectors: np.ndarray, supercell_matrix: np.ndarray, rotation_angle: float = 0.0
    ) -> np.ndarray:
        """Apply supercell matrix and rotation to lattice vectors."""
        # First apply supercell transformation
        supercell_vectors = supercell_matrix @ vectors

        # Then apply rotation if specified
        if rotation_angle != 0.0:
            rotation_matrix = self._create_rotation_matrix(rotation_angle)
            supercell_vectors = rotation_matrix @ supercell_vectors

        return supercell_vectors

    def _find_matching_substrate_supercell(
        self, film_vectors: np.ndarray, substrate_vectors: np.ndarray
    ) -> Optional[np.ndarray]:
        """Find diagonal substrate supercell matrix that matches film vectors within tolerance."""
        film_lengths = [np.linalg.norm(film_vectors[i]) for i in range(2)]
        substrate_lengths = [np.linalg.norm(substrate_vectors[i]) for i in range(2)]

        # Try diagonal matrices to match lengths
        for n in range(1, self.max_supercell_size + 1):
            for m in range(1, self.max_supercell_size + 1):
                substrate_supercell_matrix = np.array([[n, 0], [0, m]])
                substrate_supercell_vectors = substrate_supercell_matrix @ substrate_vectors
                substrate_supercell_lengths = [np.linalg.norm(substrate_supercell_vectors[i]) for i in range(2)]

                # Check if lengths match within tolerance
                length_match_a = (
                    abs(substrate_supercell_lengths[0] - film_lengths[0]) / film_lengths[0] < self.length_tolerance
                )
                length_match_b = (
                    abs(substrate_supercell_lengths[1] - film_lengths[1]) / film_lengths[1] < self.length_tolerance
                )

                if length_match_a and length_match_b:
                    # Check if vectors are parallel (aligned)
                    aligned_substrate = align_first_vector_to_x_2d_right_handed(substrate_supercell_vectors)
                    aligned_film = align_first_vector_to_x_2d_right_handed(film_vectors)

                    a_parallel = are_vectors_colinear(aligned_substrate[0], aligned_film[0])
                    b_parallel = are_vectors_colinear(aligned_substrate[1], aligned_film[1])

                    if a_parallel and b_parallel:
                        return substrate_supercell_matrix

        return None

    def _calculate_strain_matrix(self, substrate_vectors: np.ndarray, film_vectors: np.ndarray) -> np.ndarray:
        """Calculate strain matrix to transform film vectors to substrate vectors."""
        try:
            film_inv = np.linalg.inv(film_vectors)
            strain_2d = film_inv @ substrate_vectors
            return convert_2x2_to_3x3(strain_2d.tolist())
        except np.linalg.LinAlgError:
            return np.eye(3)

    def _generate_csl_matches(self) -> List[CSLMatchHolder]:
        """Generate CSL matches by testing rotations and supercell combinations."""
        if not self._validate_same_lattice_type():
            raise ValueError("Substrate and film must have the same lattice type for CSL analysis")

        substrate_vectors = self._get_lattice_vectors_2d(self.substrate_material)
        film_vectors = self._get_lattice_vectors_2d(self.film_material)

        matches = []
        match_id = 0

        # Generate supercell matrices for film
        supercell_matrices = create_2d_supercell_matrices(self.max_supercell_size)

        # Test different rotation angles
        angles = np.arange(0, self.max_rotation_angle + self.angle_step, self.angle_step)

        for film_supercell_matrix in supercell_matrices:
            film_supercell_area = abs(np.linalg.det(film_supercell_matrix))
            if film_supercell_area > self.max_area:
                continue

            for angle in angles:
                # Apply supercell and rotation to film
                film_transformed_vectors = self._apply_supercell_and_rotation(
                    film_vectors, film_supercell_matrix, angle
                )

                # Find matching substrate supercell
                substrate_supercell_matrix = self._find_matching_substrate_supercell(
                    film_transformed_vectors, substrate_vectors
                )

                if substrate_supercell_matrix is not None:
                    # Calculate substrate supercell vectors
                    substrate_supercell_vectors = substrate_supercell_matrix @ substrate_vectors

                    # Calculate strain matrix
                    strain_3d = self._calculate_strain_matrix(substrate_supercell_vectors, film_transformed_vectors)
                    strain_percentage = calculate_von_mises_strain(strain_3d)

                    # Check if strain is within tolerance
                    if strain_percentage <= self.strain_tolerance:
                        match_area = abs(np.linalg.det(substrate_supercell_matrix)) * abs(
                            np.linalg.det(substrate_vectors)
                        )

                        matches.append(
                            CSLMatchHolder(
                                match_id=match_id,
                                substrate_supercell_matrix=SupercellMatrix2DSchema(
                                    root=substrate_supercell_matrix.tolist()
                                ),
                                film_supercell_matrix=SupercellMatrix2DSchema(root=film_supercell_matrix.tolist()),
                                film_strain_matrix=Matrix3x3Schema(root=strain_3d),
                                match_area=match_area,
                                total_strain_percentage=strain_percentage,
                                film_rotation_angle=angle,
                            )
                        )
                        match_id += 1

        # Sort by strain percentage, then by area
        matches.sort(key=lambda m: (m.total_strain_percentage, m.match_area))
        return matches

    @cached_property
    def csl_match_holders(self) -> List[CSLMatchHolder]:
        """Get all CSL matches, cached for performance."""
        return self._generate_csl_matches()

    def get_strained_configuration_by_match_id(self, match_id: int) -> MatchedSubstrateFilmConfigurationHolder:
        """Get strained configuration for a specific match ID."""
        match_holders = self.csl_match_holders
        if not match_holders:
            raise ValueError("No CSL matches found. Please check the input configurations and tolerances.")
        if match_id < 0 or match_id >= len(match_holders):
            raise ValueError(f"Match ID {match_id} out of range. Available IDs: 0-{len(match_holders)-1}")

        match_holder = match_holders[match_id]
        return self._create_strained_configs_from_match(match_holder)

    def _create_strained_configs_from_match(
        self, match_holder: CSLMatchHolder
    ) -> MatchedSubstrateFilmConfigurationHolder:
        """Create strained configurations from a CSL match."""
        return self.create_matched_configuration_holder(
            substrate_slab_config=self.substrate_slab_configuration,
            film_slab_config=self.film_slab_configuration,
            match_id=match_holder.match_id,
            substrate_xy_supercell_matrix=match_holder.substrate_supercell_matrix,
            film_xy_supercell_matrix=match_holder.film_supercell_matrix,
            substrate_strain_matrix=self._no_strain_matrix,
            film_strain_matrix=match_holder.film_strain_matrix,
            total_strain_percentage=match_holder.total_strain_percentage,
        )

    def get_strained_configurations(self) -> List[MatchedSubstrateFilmConfigurationHolder]:
        """Get all strained configurations for all matches."""
        strained_configs = []

        for match_holder in self.csl_match_holders:
            config_holder = self._create_strained_configs_from_match(match_holder)
            strained_configs.append(config_holder)

        return strained_configs
