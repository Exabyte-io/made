from functools import cached_property
from typing import Optional, Tuple

import numpy as np
from mat3ra.code.entity import InMemoryEntityPydantic
from mat3ra.esse.models.core.abstract.matrix_3x3 import Matrix3x3Schema
from mat3ra.esse.models.materials_category_components.entities.auxiliary.two_dimensional.supercell_matrix_2d import (
    SupercellMatrix2DSchema,
)
from mat3ra.made.tools.build.pristine_structures.two_dimensional.slab.build_parameters import SlabBuilderParameters
from mat3ra.made.tools.build.pristine_structures.two_dimensional.slab.builder import SlabBuilder
from mat3ra.made.tools.build.pristine_structures.two_dimensional.slab.configuration import SlabConfiguration

from ...build.pristine_structures.two_dimensional.slab_strained_supercell.configuration import (
    SlabStrainedSupercellConfiguration,
)
from ...operations.core.unary import supercell
from ...utils import unwrap
from ..interface.utils.holders import MatchedSubstrateFilmConfigurationHolder
from ..utils import calculate_von_mises_strain


class InterfaceAnalyzer(InMemoryEntityPydantic):
    """
    Analyzes the interface between two slabs, calculating strain matrices and strained configurations.
    This class is used to prepare the configurations for the substrate and film materials, ensuring that the film
    material's lattice vectors match the substrate material's lattice vectors.

    `substrate_strained_configuration` and `film_strained_configuration` properties provide the strained configurations.
    """

    substrate_slab_configuration: SlabConfiguration
    film_slab_configuration: SlabConfiguration
    substrate_build_parameters: Optional[SlabBuilderParameters] = None
    film_build_parameters: Optional[SlabBuilderParameters] = None
    optimize_film_supercell: bool = False

    def get_component_material(self, configuration: SlabConfiguration):
        return SlabBuilder().get_material(configuration)

    @cached_property
    def substrate_material(self):
        return self.get_component_material(self.substrate_slab_configuration)

    @cached_property
    def film_material(self):
        return self.get_component_material(self.film_slab_configuration)

    @property
    def identity_supercell(self) -> SupercellMatrix2DSchema:
        return SupercellMatrix2DSchema(root=[[1, 0], [0, 1]])

    @property
    def _no_strain_matrix(self) -> Matrix3x3Schema:
        return Matrix3x3Schema(root=np.eye(3).tolist())

    def get_film_strain_matrix(
        self,
        substrate_lattice_vector_arrays: np.ndarray,
        film_lattice_vector_arrays: np.ndarray,
    ) -> Matrix3x3Schema:
        """Calculate film strain to match substrate using ZSL approach.

        This function calculates the 2D strain matrix needed to transform the film
        lattice vectors to match the substrate lattice vectors.

        Args:
            substrate_lattice_vector_arrays: 3x3 array of substrate lattice vectors
            film_lattice_vector_arrays: 3x3 array of film lattice vectors

        Returns:
            2x2 strain matrix that transforms film vectors to substrate vectors
        """
        substrate_vectors = np.array(substrate_lattice_vector_arrays)
        film_vectors = np.array(film_lattice_vector_arrays)

        substrate_2d_vectors = substrate_vectors[:2, :2]
        film_2d_vectors = film_vectors[:2, :2]

        try:
            inv_film = np.linalg.inv(film_2d_vectors)
        except np.linalg.LinAlgError:
            raise ValueError("Film lattice vectors are not linearly independent.")

        strain_2d = inv_film @ substrate_2d_vectors

        strain_3d = np.eye(3)
        strain_3d[:2, :2] = strain_2d
        return Matrix3x3Schema(root=strain_3d.tolist())

    def get_substrate_strain_matrix(self) -> Matrix3x3Schema:
        return self._no_strain_matrix

    def _calculate_strain_for_film_supercell(self, film_n: int, film_m: int) -> float:
        """
        Calculate strain for given film supercell configuration against unchanged substrate.

        Args:
            film_n: Film supercell multiplier in a direction
            film_m: Film supercell multiplier in b direction

        Returns:
            Von Mises strain percentage
        """

        # Apply supercell to film material
        supercell_matrix = [[film_n, 0], [0, film_m]]
        film_with_supercell = supercell(self.film_material, supercell_matrix)

        # Use existing get_film_strain_matrix with supercelled film
        strain_matrix = self.get_film_strain_matrix(
            self.substrate_material.lattice.vector_arrays, film_with_supercell.lattice.vector_arrays
        )

        return calculate_von_mises_strain(np.array(unwrap(strain_matrix.root)))

    def _find_optimal_supercell_factor_for_direction(self, substrate_length: float, film_length: float) -> int:
        """
        Finds a multiplier for the component of diagonal supercell matrix
            that will get film supercell lattice vector length around substrate length.
            The N leads to the film supercell vector length being shorter than substrate length,
            and N+1 leads to the film supercell vector length being longer than substrate length.
        """

        # Find optimal: when n*film_length < substrate_length < (n+1)*film_length
        optimal = max(1, int(substrate_length / film_length))
        if (optimal + 1) * film_length - substrate_length < substrate_length - optimal * film_length:
            optimal += 1

        return optimal

    def find_optimal_film_supercell(self) -> Tuple[int, int]:
        """
        Find optimal (n, m) supercell multipliers for film that minimize strain.

        Returns:
            Tuple of (n, m) supercell multipliers that minimize strain
        """
        substrate_vectors = np.array(self.substrate_material.lattice.vector_arrays[:2])
        film_vectors = np.array(self.film_material.lattice.vector_arrays[:2])

        substrate_lengths = [np.linalg.norm(substrate_vectors[i, :2]) for i in range(2)]
        film_lengths = [np.linalg.norm(film_vectors[i, :2]) for i in range(2)]

        optimal_values = []
        for substrate_length, film_length in zip(substrate_lengths, film_lengths):
            optimal = self._find_optimal_supercell_factor_for_direction(substrate_length, film_length)
            optimal_values.append(optimal)

        n_optimal, m_optimal = optimal_values

        # Test neighboring values to find minimum strain
        candidates = [
            (n_optimal, m_optimal),
            (n_optimal + 1, m_optimal),
            (n_optimal, m_optimal + 1),
            (n_optimal + 1, m_optimal + 1),
        ]

        min_strain = float("inf")
        best_n, best_m = n_optimal, m_optimal

        for n, m in candidates:
            strain = self._calculate_strain_for_film_supercell(n, m)
            if strain < min_strain:
                min_strain = strain
                best_n, best_m = n, m

        return best_n, best_m

    def get_component_strained_configuration(
        self,
        configuration: SlabConfiguration,
        strain_matrix: Matrix3x3Schema,
        xy_supercell_matrix: SupercellMatrix2DSchema = None,
    ) -> SlabStrainedSupercellConfiguration:
        if xy_supercell_matrix is not None:
            matrix_list = [[item.root[0], item.root[1]] for item in xy_supercell_matrix.root]
        else:
            matrix_list = [[1, 0], [0, 1]]

        return SlabStrainedSupercellConfiguration(
            stack_components=configuration.stack_components,
            direction=configuration.direction,
            strain_matrix=strain_matrix,
            xy_supercell_matrix=matrix_list,
        )

    def create_matched_configuration_holder(
        self,
        substrate_slab_config: SlabConfiguration,
        film_slab_config: SlabConfiguration,
        match_id: int = 0,
        substrate_xy_supercell_matrix: SupercellMatrix2DSchema = None,
        film_xy_supercell_matrix: SupercellMatrix2DSchema = None,
        substrate_strain_matrix: Matrix3x3Schema = None,
        film_strain_matrix: Matrix3x3Schema = None,
        total_strain_percentage: Optional[float] = None,
    ) -> MatchedSubstrateFilmConfigurationHolder:
        if substrate_strain_matrix is None:
            substrate_strain_matrix = self._no_strain_matrix
        if film_strain_matrix is None:
            film_strain_matrix = self._no_strain_matrix

        substrate_config = self.get_component_strained_configuration(
            substrate_slab_config,
            substrate_strain_matrix,
            xy_supercell_matrix=substrate_xy_supercell_matrix,
        )

        film_config = self.get_component_strained_configuration(
            film_slab_config, film_strain_matrix, xy_supercell_matrix=film_xy_supercell_matrix
        )

        return MatchedSubstrateFilmConfigurationHolder(
            match_id=match_id,
            substrate_configuration=substrate_config,
            film_configuration=film_config,
            total_strain_percentage=total_strain_percentage,
        )

    @property
    def substrate_supercell_matrix(self) -> SupercellMatrix2DSchema:
        if self.substrate_build_parameters and self.substrate_build_parameters.xy_supercell_matrix:
            return SupercellMatrix2DSchema(root=self.substrate_build_parameters.xy_supercell_matrix)
        return self.identity_supercell

    @property
    def film_supercell_matrix(self) -> SupercellMatrix2DSchema:
        if self.optimize_film_supercell:
            n, m = self.find_optimal_film_supercell()
            return SupercellMatrix2DSchema(root=[[n, 0], [0, m]])
        elif self.film_build_parameters and self.film_build_parameters.xy_supercell_matrix:
            return SupercellMatrix2DSchema(root=self.film_build_parameters.xy_supercell_matrix)
        return self.identity_supercell

    @cached_property
    def substrate_strained_configuration(self) -> SlabStrainedSupercellConfiguration:
        return self.get_component_strained_configuration(
            self.substrate_slab_configuration,
            self._no_strain_matrix,
            xy_supercell_matrix=self.substrate_supercell_matrix,
        )

    @cached_property
    def film_strained_configuration(self) -> SlabStrainedSupercellConfiguration:
        return self.get_component_strained_configuration(
            self.film_slab_configuration,
            self.get_film_strain_matrix(
                self.substrate_material.lattice.vector_arrays, self.film_material.lattice.vector_arrays
            ),
            xy_supercell_matrix=self.film_supercell_matrix,
        )
