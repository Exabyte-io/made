from functools import cached_property
from typing import List, Tuple

import numpy as np
from mat3ra.code.entity import InMemoryEntityPydantic
from mat3ra.esse.models.core.abstract.matrix_3x3 import Matrix3x3Schema
from mat3ra.esse.models.materials_category_components.entities.auxiliary.two_dimensional.supercell_matrix_2d import (
    SupercellMatrix2DSchema,
)
from pymatgen.analysis.interfaces.coherent_interfaces import ZSLGenerator, CoherentInterfaceBuilder

from mat3ra.made.tools.build.slab.builders import SlabBuilder
from mat3ra.made.tools.build.slab.configuration import (
    SlabConfiguration,
    SlabStrainedSupercellConfiguration,
)
from mat3ra.made.tools.convert import to_pymatgen


class ZSLMatchHolder(InMemoryEntityPydantic):
    film_transformation: np.ndarray
    substrate_transformation: np.ndarray


class InterfaceAnalyzer(InMemoryEntityPydantic):

    """
    Analyzes the interface between two slabs, calculating strain matrices and strained configurations.
    This class is used to prepare the configurations for the substrate and film materials, ensuring that the film
    material's lattice vectors match the substrate material's lattice vectors.

    `substrate_strained_configuration` and `film_strained_configuration` properties provide the strained configurations.
    """

    substrate_slab_configuration: SlabConfiguration
    film_slab_configuration: SlabConfiguration

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
    def identity_strain(self) -> Matrix3x3Schema:
        return Matrix3x3Schema(root=np.eye(3).tolist())

    @cached_property
    def film_strain_matrix(self) -> Matrix3x3Schema:
        substrate_vectors = np.array(self.substrate_material.lattice.vector_arrays)
        film_vectors = np.array(self.film_material.lattice.vector_arrays)

        substrate_2d_vectors = substrate_vectors[:2, :2]
        film_2d_vectors = film_vectors[:2, :2]

        try:
            inv_film = np.linalg.inv(film_2d_vectors)
        except np.linalg.LinAlgError:
            raise ValueError("Film lattice vectors are not linearly independent.")
        strain2d = inv_film @ substrate_2d_vectors

        # Convert 2D strain matrix to 3D strain matrix with identity in z-direction
        strain_3d = np.eye(3)
        strain_3d[:2, :2] = strain2d
        return Matrix3x3Schema(root=strain_3d.tolist())

    @property
    def substrate_strain_matrix(self) -> Matrix3x3Schema:
        return self.identity_strain

    def get_component_strained_configuration(
        self,
        configuration: SlabConfiguration,
        strain_matrix: Matrix3x3Schema,
    ) -> SlabStrainedSupercellConfiguration:
        return SlabStrainedSupercellConfiguration(
            stack_components=configuration.stack_components,
            direction=configuration.direction,
            strain_matrix=strain_matrix,
        )
    @property
    def substrate_supercell_matrix(self) -> SupercellMatrix2DSchema:
        return self.identity_supercell

    @property
    def film_supercell_matrix(self) -> SupercellMatrix2DSchema:
        return self.identity_supercell

    @cached_property
    def substrate_strained_configuration(self) -> SlabStrainedSupercellConfiguration:
        return self.get_component_strained_configuration(
            self.substrate_slab_configuration, self.substrate_strain_matrix
        )

    @cached_property
    def film_strained_configuration(self) -> SlabStrainedSupercellConfiguration:
        return self.get_component_strained_configuration(self.film_slab_configuration, self.film_strain_matrix)


class ZSLInterfaceAnalyzer(InterfaceAnalyzer):
    """Interface analyzer using Pymatgen's ZSL algorithm to find matching supercells."""

    max_area: float = 50.0
    max_area_ratio_tol: float = 0.09
    max_length_tol: float = 0.03
    max_angle_tol: float = 0.01

    @cached_property
    def get_pymatgen_zsl_generator(self) -> ZSLGenerator:
        return ZSLGenerator(
            max_area=self.max_area,
            max_area_ratio_tol=self.max_area_ratio_tol,
            max_length_tol=self.max_length_tol,
            max_angle_tol=self.max_angle_tol,
        )

    @cached_property
    def get_pymatgen_coherent_interface_builder(self) -> CoherentInterfaceBuilder:
        return CoherentInterfaceBuilder(
            substrate_structure=to_pymatgen(self.substrate_slab_configuration.atomic_layers.crystal),
            film_structure=to_pymatgen(self.film_slab_configuration.atomic_layers.crystal),
            substrate_miller=self.substrate_slab_configuration.atomic_layers.miller_indices,
            film_miller=self.film_slab_configuration.atomic_layers.miller_indices,
            zslgen=self.get_pymatgen_zsl_generator,
        )

    @cached_property
    def zsl_match_holders(self) -> List[ZSLMatchHolder]:
        """Get ZSL matches between substrate and film slabs."""
        zsl_matches = self.get_pymatgen_coherent_interface_builder.zsl_matches
        match_holders = []
        for match in zsl_matches:
            match_holder = ZSLMatchHolder(
                film_transformation=match.film_transformation,
                substrate_transformation=match.substrate_transformation,
            )
            match_holders.append(match_holder)
        return match_holders

    def _create_strained_configs(
        self, match: ZSLMatchHolder
    ) -> Tuple[SlabStrainedSupercellConfiguration, SlabStrainedSupercellConfiguration]:
        # Use ZSL transformation matrices as supercell matrices
        substrate_supercell_matrix = match.substrate_transformation.astype(int)
        film_supercell_matrix = match.film_transformation.astype(int)

        # Substrate remains unstrained (identity strain)
        substrate_strain_2d = np.eye(2)

        # Calculate film strain using our existing approach (from InterfaceAnalyzer)
        film_strain_2d = self._calculate_film_strain_for_supercells(substrate_supercell_matrix, film_supercell_matrix)

        substrate_supercell_2d = SupercellMatrix2DSchema(root=substrate_supercell_matrix.tolist())
        film_supercell_2d = SupercellMatrix2DSchema(root=film_supercell_matrix.tolist())

        # Convert 2D strain matrices to 3D
        substrate_strain_3d = np.eye(3)
        substrate_strain_3d[:2, :2] = substrate_strain_2d

        film_strain_3d = np.eye(3)
        film_strain_3d[:2, :2] = film_strain_2d

        substrate_config = SlabStrainedSupercellConfiguration(
            stack_components=self.substrate_slab_configuration.stack_components,
            direction=self.substrate_slab_configuration.direction,
            xy_supercell_matrix=substrate_supercell_2d,
            strain_matrix=Matrix3x3Schema(root=substrate_strain_3d.tolist()),
        )

        film_config = SlabStrainedSupercellConfiguration(
            stack_components=self.film_slab_configuration.stack_components,
            direction=self.film_slab_configuration.direction,
            xy_supercell_matrix=film_supercell_2d,
            strain_matrix=Matrix3x3Schema(root=film_strain_3d.tolist()),
        )

        return substrate_config, film_config

    def _calculate_film_strain_for_supercells(
        self, substrate_supercell_matrix: np.ndarray, film_supercell_matrix: np.ndarray
    ) -> np.ndarray:
        """Calculate film strain to match substrate using our existing approach."""
        # Get original 2D lattice vectors
        substrate_original_2d = np.array(self.substrate_material.lattice.vector_arrays[:2])[:, :2]
        film_original_2d = np.array(self.film_material.lattice.vector_arrays[:2])[:, :2]

        # Calculate supercell lattices
        substrate_supercell_2d = substrate_original_2d @ substrate_supercell_matrix
        film_supercell_2d = film_original_2d @ film_supercell_matrix

        # Calculate strain to match substrate supercell
        try:
            inv_film_supercell = np.linalg.inv(film_supercell_2d)
        except np.linalg.LinAlgError:
            raise ValueError("Film supercell lattice vectors are not linearly independent.")

        strain_2d = inv_film_supercell @ substrate_supercell_2d
        return strain_2d

    def get_strained_configurations(
        self,
    ) -> List[Tuple[SlabStrainedSupercellConfiguration, SlabStrainedSupercellConfiguration]]:
        """Get strained configurations for all ZSL matches."""
        strained_configs = []
        for match in self.zsl_match_holders:
            substrate_config, film_config = self._create_strained_configs(match)
            strained_configs.append((substrate_config, film_config))
        return strained_configs
