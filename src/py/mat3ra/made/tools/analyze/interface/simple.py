from functools import cached_property
from typing import Optional

import numpy as np
from mat3ra.code.entity import InMemoryEntityPydantic
from mat3ra.esse.models.core.abstract.matrix_3x3 import Matrix3x3Schema
from mat3ra.esse.models.materials_category_components.entities.auxiliary.two_dimensional.supercell_matrix_2d import (
    SupercellMatrix2DSchema,
)

from ...build.slab.builders import SlabBuilder, SlabBuilderParameters
from ...build.slab.configurations import SlabConfiguration, SlabStrainedSupercellConfiguration
from ..interface.utils.holders import MatchedSubstrateFilmConfigurationHolder


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
        if self.film_build_parameters and self.film_build_parameters.xy_supercell_matrix:
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
