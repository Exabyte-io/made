from functools import cached_property
from functools import cached_property
from typing import List, Dict

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


class ZSLConfigurationWithStrain(InMemoryEntityPydantic):
    """Configuration with ZSL strain metadata."""

    substrate_config: SlabStrainedSupercellConfiguration
    film_config: SlabStrainedSupercellConfiguration
    strain_info: Dict[str, float]


class ZSLInterfaceAnalyzer(BaseInterfaceAnalyzer):
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
        zsl_matches = self.get_pymatgen_coherent_interface_builder.zsl_matches()
        match_holders = []
        for match in zsl_matches:
            film_transformation = match.film_transformation
            substrate_transformation = match.substrate_transformation
            film_sl_vectors = match.film_sl_vectors
            substrate_sl_vectors = match.substrate_sl_vectors

            match_holder = ZSLMatchHolder(
                film_transformation=film_transformation,
                substrate_transformation=substrate_transformation,
                film_sl_vectors=film_sl_vectors,
                substrate_sl_vectors=substrate_sl_vectors,
            )
            match_holders.append(match_holder)
        return match_holders

    def _create_strained_config(self, match: ZSLMatchHolder) -> SlabStrainedSupercellConfiguration:
        # break transformations into 2D supercell and strain matrices

