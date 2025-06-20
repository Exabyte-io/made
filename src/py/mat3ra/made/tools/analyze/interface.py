from functools import cached_property
from typing import List, Tuple, Dict, Any

import numpy as np
from mat3ra.code.entity import InMemoryEntityPydantic
from mat3ra.esse.models.core.abstract.matrix_3x3 import Matrix3x3Schema
from mat3ra.esse.models.materials_category_components.entities.auxiliary.two_dimensional.supercell_matrix_2d import (
    SupercellMatrix2DSchema,
)
from pymatgen.analysis.interfaces.coherent_interfaces import (
    ZSLGenerator,
)

from mat3ra.made.tools.build.slab.builders import SlabBuilder
from mat3ra.made.tools.build.slab.configuration import SlabConfiguration, SlabStrainedSupercellConfiguration
from mat3ra.made.tools.build.interface.enums import StrainModes
from ..convert import to_pymatgen
from ..modify import translate_to_z_level


class InterfaceAnalyzer(InMemoryEntityPydantic):
    substrate_slab_configuration: SlabConfiguration
    film_slab_configuration: SlabConfiguration

    @cached_property
    def substrate_material(self):
        return SlabBuilder().get_material(self.substrate_slab_configuration)

    @cached_property
    def film_material(self):
        return SlabBuilder().get_material(self.film_slab_configuration)

    @cached_property
    def film_strain_matrix(self) -> Matrix3x3Schema:
        substrate_vectors = np.array(self.substrate_material.lattice.vector_arrays)
        film_vectors = np.array(self.film_material.lattice.vector_arrays)

        substrate_2d_vectors = substrate_vectors[:2, :2]
        film_2d_vectors = film_vectors[:2, :2]

        try:
            inv_film_2d_vectors = np.linalg.inv(film_2d_vectors)
        except np.linalg.LinAlgError:
            raise ValueError("Film lattice vectors are not linearly independent.")

        strain_2d_matrix = inv_film_2d_vectors @ substrate_2d_vectors

        film_strain_matrix_3x3 = np.identity(3)
        film_strain_matrix_3x3[:2, :2] = strain_2d_matrix
        return Matrix3x3Schema(root=film_strain_matrix_3x3.tolist())

    @property
    def substrate_strain_matrix(self) -> Matrix3x3Schema:
        return Matrix3x3Schema(root=np.identity(3).tolist())

    @property
    def substrate_supercell_matrix(self) -> SupercellMatrix2DSchema:
        return SupercellMatrix2DSchema(root=[[1.0, 0.0], [0.0, 1.0]])

    @property
    def film_supercell_matrix(self) -> SupercellMatrix2DSchema:
        return SupercellMatrix2DSchema(root=[[1.0, 0.0], [0.0, 1.0]])

    @cached_property
    def substrate_strained_configuration(self) -> SlabStrainedSupercellConfiguration:
        return SlabStrainedSupercellConfiguration(
            stack_components=self.substrate_slab_configuration.stack_components,
            direction=self.substrate_slab_configuration.direction,
            xy_supercell_matrix=self.substrate_supercell_matrix,
            strain_matrix=self.substrate_strain_matrix,
        )

    @cached_property
    def film_strained_configuration(self) -> SlabStrainedSupercellConfiguration:
        return SlabStrainedSupercellConfiguration(
            stack_components=self.film_slab_configuration.stack_components,
            direction=self.film_slab_configuration.direction,
            xy_supercell_matrix=self.film_supercell_matrix,
            strain_matrix=self.film_strain_matrix,
        )


class ZSLConfigurationWithStrain(InMemoryEntityPydantic):
    """Configuration with ZSL strain information."""

    substrate_config: SlabStrainedSupercellConfiguration
    film_config: SlabStrainedSupercellConfiguration
    strain_info: Dict[str, float]


class ZSLInterfaceAnalyzer(InMemoryEntityPydantic):
    substrate_slab_configuration: SlabConfiguration
    film_slab_configuration: SlabConfiguration
    max_area: float = 50.0
    max_area_ratio_tol: float = 0.09
    max_length_tol: float = 0.03
    max_angle_tol: float = 0.01

    @cached_property
    def substrate_material(self):
        return SlabBuilder().get_material(self.substrate_slab_configuration)

    @cached_property
    def film_material(self):
        return SlabBuilder().get_material(self.film_slab_configuration)

    @cached_property
    def _zsl_generator(self) -> ZSLGenerator:
        return ZSLGenerator(
            max_area=self.max_area,
            max_area_ratio_tol=self.max_area_ratio_tol,
            max_length_tol=self.max_length_tol,
            max_angle_tol=self.max_angle_tol,
        )

    def _extract_strain_from_zsl_match(self, match) -> Dict[str, float]:
        """Extract strain information from a ZSL match."""
        # Get the original 2D lattice vectors (in-plane vectors a and b)
        substrate_vectors = np.array(self.substrate_material.lattice.vector_arrays)
        film_vectors = np.array(self.film_material.lattice.vector_arrays)
        substrate_2d = substrate_vectors[:2, :2]  # Take first 2 rows and first 2 columns
        film_2d = film_vectors[:2, :2]  # Take first 2 rows and first 2 columns

        # Get the superlattice vectors from the ZSL match
        # These are the matched superlattice vectors that should have nearly the same lengths and angles
        substrate_sl_2d = np.array(match.substrate_sl_vectors)  # These are the actual superlattice vectors
        film_sl_2d = np.array(match.film_sl_vectors)  # These are the actual superlattice vectors

        # Calculate strains by comparing superlattice vectors to original vectors
        # For substrate strain: (|sl_vector| - |original_vector|) / |original_vector|
        substrate_a_strain = (np.linalg.norm(substrate_sl_2d[0]) - np.linalg.norm(substrate_2d[0])) / np.linalg.norm(
            substrate_2d[0]
        )
        substrate_b_strain = (np.linalg.norm(substrate_sl_2d[1]) - np.linalg.norm(substrate_2d[1])) / np.linalg.norm(
            substrate_2d[1]
        )

        film_a_strain = (np.linalg.norm(film_sl_2d[0]) - np.linalg.norm(film_2d[0])) / np.linalg.norm(film_2d[0])
        film_b_strain = (np.linalg.norm(film_sl_2d[1]) - np.linalg.norm(film_2d[1])) / np.linalg.norm(film_2d[1])

        # Calculate average strains
        substrate_avg_strain = (abs(substrate_a_strain) + abs(substrate_b_strain)) / 2
        film_avg_strain = (abs(film_a_strain) + abs(film_b_strain)) / 2

        return {
            "substrate_strain": substrate_avg_strain,
            "film_strain": film_avg_strain,
            "substrate_a_strain": substrate_a_strain,
            "substrate_b_strain": substrate_b_strain,
            "film_a_strain": film_a_strain,
            "film_b_strain": film_b_strain,
        }

    def _get_supercell_matrix_from_sl_vectors(self, sl_vectors) -> np.ndarray:
        """
        Convert superlattice vectors to supercell matrix.

        Args:
            sl_vectors: List of superlattice vectors from ZSL match

        Returns:
            2x2 supercell matrix as numpy array
        """
        # The superlattice vectors are already the transformation matrix
        # Convert to numpy array and ensure it's 2x2
        matrix = np.array(sl_vectors)
        if matrix.shape != (2, 2):
            raise ValueError(f"Expected 2x2 superlattice vectors, got shape {matrix.shape}")

        # Round to nearest integers as supercell matrices should be integers
        return np.round(matrix).astype(int)

    def get_strained_slab_configurations(
        self,
    ) -> List[Tuple[SlabStrainedSupercellConfiguration, SlabStrainedSupercellConfiguration]]:
        """
        Generate strained slab configurations for substrate and film using ZSL algorithm.

        Returns:
            List of tuples (substrate_config, film_config) for each ZSL match.
        """
        from pymatgen.analysis.interfaces.zsl import vec_area

        substrate_vectors = np.array(self.substrate_material.lattice.vector_arrays)
        film_vectors = np.array(self.film_material.lattice.vector_arrays)
        substrate_2d = substrate_vectors[:2, :2]
        film_2d = film_vectors[:2, :2]

        substrate_area = vec_area(substrate_2d[0], substrate_2d[1])
        film_area = vec_area(film_2d[0], film_2d[1])

        configurations = []

        zsl_generator = self._zsl_generator

        transformation_sets = zsl_generator.generate_sl_transformation_sets(film_area, substrate_area)

        zsl_matches = zsl_generator.get_equiv_transformations(transformation_sets, film_2d, substrate_2d)

        for match in zsl_matches:
            # match is a list with 4 arrays: [film_sl_vectors, substrate_sl_vectors, film_transformation, substrate_transformation]
            film_sl_vectors = match[0]
            substrate_sl_vectors = match[1]
            film_transformation = match[2]
            substrate_transformation = match[3]

            substrate_transformed_2d = substrate_2d @ substrate_transformation
            try:
                substrate_strain_2d = substrate_sl_vectors @ np.linalg.inv(substrate_transformed_2d)
            except np.linalg.LinAlgError:
                substrate_strain_2d = np.eye(2)

            film_transformed_2d = film_2d @ film_transformation
            try:
                film_strain_2d = substrate_sl_vectors @ np.linalg.inv(film_transformed_2d)
            except np.linalg.LinAlgError:
                film_strain_2d = np.eye(2)

            # Convert 2D strain matrices to 3D
            substrate_strain_3d = np.eye(3)
            substrate_strain_3d[:2, :2] = substrate_strain_2d

            film_strain_3d = np.eye(3)
            film_strain_3d[:2, :2] = film_strain_2d

            # Create strained slab configurations
            substrate_config = SlabStrainedSupercellConfiguration(
                stack_components=self.substrate_slab_configuration.stack_components,
                direction=self.substrate_slab_configuration.direction,
                xy_supercell_matrix=substrate_transformation.astype(int).tolist(),
                strain_matrix=substrate_strain_3d.tolist(),
            )
            film_config = SlabStrainedSupercellConfiguration(
                stack_components=self.film_slab_configuration.stack_components,
                direction=self.film_slab_configuration.direction,
                xy_supercell_matrix=film_transformation.astype(int).tolist(),
                strain_matrix=film_strain_3d.tolist(),
            )

            configurations.append((substrate_config, film_config))

        return configurations

    def get_strained_slab_configurations_with_metadata(self) -> List[ZSLConfigurationWithStrain]:
        """
        Generate strained slab configurations with strain metadata.

        Returns:
            List of ZSLConfigurationWithStrain objects containing configurations and strain info.
        """
        from pymatgen.analysis.interfaces.zsl import vec_area

        # Get the 2D (in-plane) lattice vectors for ZSL analysis
        substrate_vectors = np.array(self.substrate_material.lattice.vector_arrays)
        film_vectors = np.array(self.film_material.lattice.vector_arrays)
        substrate_2d = substrate_vectors[:2, :2]  # Take first 2 rows and first 2 columns
        film_2d = film_vectors[:2, :2]  # Take first 2 rows and first 2 columns

        # Calculate areas of the 2D lattices
        substrate_area = vec_area(substrate_2d[0], substrate_2d[1])
        film_area = vec_area(film_2d[0], film_2d[1])

        configurations = []

        # Use the cached ZSL generator
        zsl_generator = self._zsl_generator

        # Generate transformation sets for the given areas
        transformation_sets = zsl_generator.generate_sl_transformation_sets(film_area, substrate_area)

        # Get equivalent transformations that produce matching superlattices
        zsl_matches = zsl_generator.get_equiv_transformations(transformation_sets, film_2d, substrate_2d)

        for match in zsl_matches:
            # match is a list with 4 arrays: [film_sl_vectors, substrate_sl_vectors, film_transformation, substrate_transformation]
            film_sl_vectors = match[0]
            substrate_sl_vectors = match[1]
            film_transformation = match[2]
            substrate_transformation = match[3]

            # Calculate strain matrices to transform original lattices to matching superlattices
            substrate_transformed_2d = substrate_2d @ substrate_transformation
            try:
                substrate_strain_2d = substrate_sl_vectors @ np.linalg.inv(substrate_transformed_2d)
            except np.linalg.LinAlgError:
                substrate_strain_2d = np.eye(2)

            film_transformed_2d = film_2d @ film_transformation
            try:
                film_strain_2d = substrate_sl_vectors @ np.linalg.inv(film_transformed_2d)
            except np.linalg.LinAlgError:
                film_strain_2d = np.eye(2)

            # Convert 2D strain matrices to 3D
            substrate_strain_3d = np.eye(3)
            substrate_strain_3d[:2, :2] = substrate_strain_2d

            film_strain_3d = np.eye(3)
            film_strain_3d[:2, :2] = film_strain_2d

            # Create strained slab configurations
            substrate_config = SlabStrainedSupercellConfiguration(
                stack_components=self.substrate_slab_configuration.stack_components,
                direction=self.substrate_slab_configuration.direction,
                xy_supercell_matrix=substrate_transformation.astype(int).tolist(),
                strain_matrix=substrate_strain_3d.tolist(),
            )
            film_config = SlabStrainedSupercellConfiguration(
                stack_components=self.film_slab_configuration.stack_components,
                direction=self.film_slab_configuration.direction,
                xy_supercell_matrix=film_transformation.astype(int).tolist(),
                strain_matrix=film_strain_3d.tolist(),
            )

            # Calculate strain information using the strain matrices
            strain_info = self._extract_strain_from_strain_matrices(substrate_strain_2d, film_strain_2d)

            configurations.append(
                ZSLConfigurationWithStrain(
                    substrate_config=substrate_config, film_config=film_config, strain_info=strain_info
                )
            )

        # Sort by average strain (lowest first)
        configurations.sort(key=lambda x: (x.strain_info["substrate_strain"] + x.strain_info["film_strain"]) / 2)

        return configurations

    def _extract_strain_from_strain_matrices(self, substrate_strain_2d, film_strain_2d) -> Dict[str, float]:
        """Extract strain information from strain matrices."""
        # Calculate strains from strain matrices
        # The strain matrix represents the deformation: strained = strain_matrix @ original
        # For strain calculation: strain = (strained_length - original_length) / original_length
        # If strain_matrix is S, then strain = |S| - 1 for each principal direction

        # Calculate the determinant and eigenvalues to get strain information
        substrate_det = np.linalg.det(substrate_strain_2d)
        film_det = np.linalg.det(film_strain_2d)

        # Get eigenvalues to understand principal strains
        substrate_eigenvals = np.linalg.eigvals(substrate_strain_2d)
        film_eigenvals = np.linalg.eigvals(film_strain_2d)

        # Calculate strains as (stretch_factor - 1)
        substrate_strains = np.abs(substrate_eigenvals - 1.0)
        film_strains = np.abs(film_eigenvals - 1.0)

        # Average strain for each material
        substrate_avg_strain = np.mean(substrate_strains)
        film_avg_strain = np.mean(film_strains)

        return {
            "substrate_strain": float(substrate_avg_strain),
            "film_strain": float(film_avg_strain),
            "substrate_a_strain": float(substrate_strains[0]),
            "substrate_b_strain": float(substrate_strains[1]),
            "film_a_strain": float(film_strains[0]),
            "film_b_strain": float(film_strains[1]),
        }
