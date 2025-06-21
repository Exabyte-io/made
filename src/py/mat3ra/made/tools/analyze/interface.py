from functools import cached_property
from functools import cached_property
from typing import List, Dict, Tuple

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
    film_sl_vectors: np.ndarray
    substrate_sl_vectors: np.ndarray
    match_transformation: np.ndarray


class BaseInterfaceAnalyzer(InMemoryEntityPydantic):
    substrate_slab_configuration: SlabConfiguration
    film_slab_configuration: SlabConfiguration

    @cached_property
    def substrate_material(self):
        return SlabBuilder().get_material(self.substrate_slab_configuration)

    @cached_property
    def film_material(self):
        return SlabBuilder().get_material(self.film_slab_configuration)

    @property
    def identity_supercell(self) -> SupercellMatrix2DSchema:
        return SupercellMatrix2DSchema(root=[[1, 0], [0, 1]])

    @property
    def identity_strain(self) -> Matrix3x3Schema:
        return Matrix3x3Schema(root=np.eye(3).tolist())

    def to_2d(self, lattice) -> np.ndarray:
        """Extract 2D lattice vectors (a and b vectors) from 3D lattice."""
        return np.array(lattice.vector_arrays[:2])[:, :2]

    def to_3d_strain(self, strain_2d: np.ndarray) -> Matrix3x3Schema:
        """Convert 2D strain matrix to 3D strain matrix with identity in z-direction."""
        strain_3d = np.eye(3)
        strain_3d[:2, :2] = strain_2d
        return Matrix3x3Schema(root=strain_3d.tolist())


class InterfaceAnalyzer(BaseInterfaceAnalyzer):
    """Simple interface analyzer without ZSL: computes direct strain and identity supercell."""

    @cached_property
    def film_strain_matrix(self) -> Matrix3x3Schema:
        sub2d = self.to_2d(self.substrate_material.lattice)
        film2d = self.to_2d(self.film_material.lattice)
        try:
            inv_film = np.linalg.inv(film2d)
        except np.linalg.LinAlgError:
            raise ValueError("Film lattice vectors are not linearly independent.")
        strain2d = inv_film @ sub2d
        return self.to_3d_strain(strain2d)

    @property
    def substrate_strain_matrix(self) -> Matrix3x3Schema:
        return self.identity_strain

    @property
    def substrate_supercell_matrix(self) -> SupercellMatrix2DSchema:
        return self.identity_supercell

    @property
    def film_supercell_matrix(self) -> SupercellMatrix2DSchema:
        return self.identity_supercell

    @cached_property
    def substrate_strained_configuration(self) -> SlabStrainedSupercellConfiguration:
        return SlabStrainedSupercellConfiguration(
            stack_components=self.substrate_slab_configuration.stack_components,
            direction=self.substrate_slab_configuration.direction,
            xy_supercell_matrix=self.identity_supercell,
            strain_matrix=self.identity_strain,
        )

    @cached_property
    def film_strained_configuration(self) -> SlabStrainedSupercellConfiguration:
        return SlabStrainedSupercellConfiguration(
            stack_components=self.film_slab_configuration.stack_components,
            direction=self.film_slab_configuration.direction,
            xy_supercell_matrix=self.identity_supercell,
            strain_matrix=self.film_strain_matrix,
        )


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
            film_transformation = match.film_transformation
            substrate_transformation = match.substrate_transformation
            film_sl_vectors = match.film_sl_vectors
            substrate_sl_vectors = match.substrate_sl_vectors
            match_transformation = match.match_transformation

            match_holder = ZSLMatchHolder(
                film_transformation=film_transformation,
                substrate_transformation=substrate_transformation,
                film_sl_vectors=film_sl_vectors,
                substrate_sl_vectors=substrate_sl_vectors,
                match_transformation=match_transformation,
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
        film_strain_2d = self._calculate_film_strain_for_supercells(
            substrate_supercell_matrix, film_supercell_matrix
        )
        
        substrate_supercell_2d = SupercellMatrix2DSchema(root=substrate_supercell_matrix.tolist())
        film_supercell_2d = SupercellMatrix2DSchema(root=film_supercell_matrix.tolist())
        
        substrate_strain_3d = self.to_3d_strain(substrate_strain_2d)
        film_strain_3d = self.to_3d_strain(film_strain_2d)

        substrate_config = SlabStrainedSupercellConfiguration(
            stack_components=self.substrate_slab_configuration.stack_components,
            direction=self.substrate_slab_configuration.direction,
            xy_supercell_matrix=substrate_supercell_2d,
            strain_matrix=substrate_strain_3d,
        )

        film_config = SlabStrainedSupercellConfiguration(
            stack_components=self.film_slab_configuration.stack_components,
            direction=self.film_slab_configuration.direction,
            xy_supercell_matrix=film_supercell_2d,
            strain_matrix=film_strain_3d,
        )

        return substrate_config, film_config

    def _calculate_film_strain_for_supercells(
        self, substrate_supercell_matrix: np.ndarray, film_supercell_matrix: np.ndarray
    ) -> np.ndarray:
        """Calculate film strain to match substrate using our existing approach."""
        # Get original 2D lattice vectors
        substrate_original_2d = self.to_2d(self.substrate_material.lattice)
        film_original_2d = self.to_2d(self.film_material.lattice)
        
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
