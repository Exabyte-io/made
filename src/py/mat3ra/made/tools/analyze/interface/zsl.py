from typing import List

import numpy as np
from mat3ra.code.entity import InMemoryEntityPydantic
from mat3ra.esse.models.core.abstract.matrix_3x3 import Matrix3x3Schema
from mat3ra.esse.models.materials_category_components.entities.auxiliary.two_dimensional.supercell_matrix_2d import (
    SupercellMatrix2DSchema,
)
from mat3ra.utils.matrix import convert_2x2_to_3x3
from pymatgen.analysis.interfaces.coherent_interfaces import (
    ZSLGenerator as PymatgenZSLGenerator,
)

from mat3ra.made.tools.analyze.interface.simple import InterfaceAnalyzer
from mat3ra.made.tools.analyze.interface.utils.holders import MatchedSubstrateFilmConfigurationHolder
from mat3ra.made.tools.build import MaterialWithBuildMetadata
from mat3ra.made.tools.build.slab.builders import SlabBuilder
from mat3ra.made.utils import calculate_von_mises_strain


class ZSLMatchHolder(InMemoryEntityPydantic):
    match_id: int
    substrate_supercell_matrix: SupercellMatrix2DSchema
    film_supercell_matrix: SupercellMatrix2DSchema
    strain_transformation_matrix: Matrix3x3Schema
    substrate_sl_vector_arrays: np.ndarray = None
    film_sl_vector_arrays: np.ndarray = None
    film_vectors: np.ndarray = None
    substrate_vectors: np.ndarray = None
    match_area: float
    total_strain_percentage: float


class ZSLInterfaceAnalyzer(InterfaceAnalyzer):
    """Interface analyzer using Pymatgen's ZSL algorithm to find matching supercells."""

    max_area: float = 50.0
    max_area_ratio_tol: float = 0.09
    max_length_tol: float = 0.03
    max_angle_tol: float = 0.01

    @classmethod
    def calculate_total_strain_percentage(cls, strain_matrix: list) -> float:
        return calculate_von_mises_strain(strain_matrix)

    @property
    def film_slab(self) -> MaterialWithBuildMetadata:
        return SlabBuilder().get_material(self.film_slab_configuration)

    @property
    def substrate_slab(self) -> MaterialWithBuildMetadata:
        return SlabBuilder().get_material(self.substrate_slab_configuration)

    def get_pymatgen_match_holders(self):
        zsl_generator = PymatgenZSLGenerator(
            max_area=self.max_area,
            max_area_ratio_tol=self.max_area_ratio_tol,
            max_length_tol=self.max_length_tol,
            max_angle_tol=self.max_angle_tol,
        )

        film_vectors = self.film_slab.lattice.vector_arrays[0:2]
        substrate_vectors = self.substrate_slab.lattice.vector_arrays[0:2]

        zsl_matches_pymatgen = list(
            zsl_generator(
                film_vectors=film_vectors,
                substrate_vectors=substrate_vectors,
            )
        )
        return zsl_matches_pymatgen

    @property
    def zsl_match_holders(self) -> List[ZSLMatchHolder]:
        match_holders = []
        for idx, match_pymatgen in enumerate(self.get_pymatgen_match_holders()):
            real_film_supercell_matrix = (
                match_pymatgen.film_sl_vectors[:, :2] @ np.linalg.inv(np.array(match_pymatgen.film_vectors)[:, :2])
            ).astype(int)

            real_substrate_supercell_matrix = (
                match_pymatgen.substrate_sl_vectors[:, :2]
                @ np.linalg.inv(np.array(match_pymatgen.substrate_vectors)[:, :2])
            ).astype(int)

            film_vectors = np.array(self.film_slab.lattice.vector_arrays[0:2])[:, :2]
            substrate_vectors = np.array(self.substrate_slab.lattice.vector_arrays[0:2])[:, :2]

            film_sl_supercell_vectors = film_vectors @ real_film_supercell_matrix
            substrate_sl_supercell_vectors = substrate_vectors @ real_substrate_supercell_matrix
            if (
                abs(np.linalg.det(film_sl_supercell_vectors)) < 1e-4
                or abs(np.linalg.det(substrate_sl_supercell_vectors)) < 1e-4
            ):
                continue
            real_strain_matrix = np.linalg.inv(film_sl_supercell_vectors) @ substrate_sl_supercell_vectors

            # Sanity check for strain matrix: applied to film sl vectors should yield substrate sl vectors
            delta = np.abs(
                film_vectors @ real_film_supercell_matrix @ real_strain_matrix
                - substrate_vectors @ real_substrate_supercell_matrix
            )
            assert np.allclose(delta, 0.0, atol=1e-5)

            real_strain_matrix = convert_2x2_to_3x3(real_strain_matrix)

            match_holder = ZSLMatchHolder(
                match_id=idx,
                substrate_supercell_matrix=SupercellMatrix2DSchema(root=real_substrate_supercell_matrix.tolist()),
                film_supercell_matrix=SupercellMatrix2DSchema(root=real_film_supercell_matrix.tolist()),
                strain_transformation_matrix=Matrix3x3Schema(root=real_strain_matrix),
                match_area=match_pymatgen.match_area,
                total_strain_percentage=self.calculate_total_strain_percentage(real_strain_matrix),
            )
            match_holders.append(match_holder)

        # sort matches by area in ascending order
        match_holders.sort(key=lambda x: x.match_area)
        # sort matches by strain in ascending order
        match_holders.sort(key=lambda x: x.total_strain_percentage)

        return match_holders

    def get_strained_configuration_by_match_id(self, match_id: int) -> MatchedSubstrateFilmConfigurationHolder:
        match_holders = self.zsl_match_holders
        if not match_holders:
            raise ValueError("No ZSL matches found. Please check the input configurations.")
        if match_id < 0 or match_id >= len(match_holders):
            raise ValueError(f"Match ID {match_id} out of range. Available IDs: 0-{len(match_holders)-1}")

        match_holder = match_holders[match_id]
        return self._create_strained_configs_from_match(match_holder)

    def _create_strained_configs_from_match(
        self, match_holder: ZSLMatchHolder
    ) -> MatchedSubstrateFilmConfigurationHolder:

        return super().create_matched_configuration_holder(
            substrate_slab_config=self.substrate_slab_configuration,
            film_slab_config=self.film_slab_configuration,
            match_id=match_holder.match_id,
            substrate_xy_supercell_matrix=match_holder.substrate_supercell_matrix,
            film_xy_supercell_matrix=match_holder.film_supercell_matrix,
            substrate_strain_matrix=self._no_strain_matrix,
            film_strain_matrix=Matrix3x3Schema(root=match_holder.strain_transformation_matrix.root),
            total_strain_percentage=match_holder.total_strain_percentage,
        )

    def get_strained_configurations(
        self,
    ) -> List[MatchedSubstrateFilmConfigurationHolder]:
        strained_configs = []

        for match_holder in self.zsl_match_holders:
            config_holder = self._create_strained_configs_from_match(match_holder)
            strained_configs.append(config_holder)

        return strained_configs
