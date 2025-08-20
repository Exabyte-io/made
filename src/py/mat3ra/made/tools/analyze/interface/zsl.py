from collections import defaultdict
from functools import cached_property
from typing import ClassVar, List, Optional

import numpy as np
from mat3ra.code.entity import InMemoryEntityPydantic
from mat3ra.esse.models.core.abstract.matrix_3x3 import Matrix3x3Schema
from mat3ra.esse.models.materials_category_components.entities.auxiliary.two_dimensional.supercell_matrix_2d import (
    SupercellMatrix2DSchema,
)
from mat3ra.made.tools.build.pristine_structures.two_dimensional.slab.builder import SlabBuilder
from mat3ra.utils.matrix import convert_2x2_to_3x3
from pymatgen.analysis.interfaces.coherent_interfaces import ZSLGenerator as PymatgenZSLGenerator

from ...build_components import MaterialWithBuildMetadata
from ..interface.simple import InterfaceAnalyzer
from ..interface.utils.holders import MatchedSubstrateFilmConfigurationHolder
from ..utils import calculate_von_mises_strain
from .utils.vector import get_global_gcd


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

    max_area: float = 400.0
    max_area_ratio_tol: float = 0.09
    max_length_tol: float = 0.03
    max_angle_tol: float = 0.01
    reduce_result_cell: bool = True
    bidirectional: bool = True

    math_precision: ClassVar[float] = 1e-4

    @classmethod
    def calculate_total_strain_percentage(cls, strain_matrix: list) -> float:
        return calculate_von_mises_strain(np.array(strain_matrix))

    @cached_property
    def film_slab(self) -> MaterialWithBuildMetadata:
        return SlabBuilder().get_material(self.film_slab_configuration)

    @cached_property
    def substrate_slab(self) -> MaterialWithBuildMetadata:
        return SlabBuilder().get_material(self.substrate_slab_configuration)

    @cached_property
    def generated_matches(self):
        # Be tolerant to pmg versions that don't support 'bidirectional'
        kwargs = dict(
            max_area=self.max_area,
            max_area_ratio_tol=self.max_area_ratio_tol,
            max_length_tol=self.max_length_tol,
            max_angle_tol=self.max_angle_tol,
        )
        try:
            # Prefer bidirectional search when available (finds more “simple” matches)
            kwargs["bidirectional"] = getattr(self, "bidirectional", True)
            zsl_generator = PymatgenZSLGenerator(**kwargs)
        except TypeError:
            # Older pymatgen: no 'bidirectional' kw
            kwargs.pop("bidirectional", None)
            zsl_generator = PymatgenZSLGenerator(**kwargs)

        film_vectors = self.film_slab.lattice.vector_arrays[0:2]
        substrate_vectors = self.substrate_slab.lattice.vector_arrays[0:2]

        return list(
            zsl_generator(
                film_vectors=film_vectors,
                substrate_vectors=substrate_vectors,
            )
        )

    @cached_property
    def zsl_match_holders(self) -> List[ZSLMatchHolder]:
        match_holders = []
        for idx, match_pymatgen in enumerate(self.generated_matches):
            match_holder = self.convert_generated_match_to_match_holder(idx, match_pymatgen)
            if match_holder is None:
                continue
            match_holders.append(match_holder)

        # sort matches by strain in ascending order, then for each equal strain by area in ascending order
        match_holders = self.sort_by_strain_then_area(match_holders)

        return match_holders

    def convert_generated_match_to_match_holder(self, match_id: int, match_pymatgen) -> Optional[ZSLMatchHolder]:
        # In-plane bases in original frames (no ad-hoc alignment)
        film_slab_vectors = np.array(self.film_slab.lattice.vector_arrays[0:2])[:, :2]  # (2,2)
        substrate_slab_vectors = np.array(self.substrate_slab.lattice.vector_arrays[0:2])[:, :2]  # (2,2)

        # Superlattice vectors from pymatgen ZSLMatch
        pmg_film_SL = match_pymatgen.film_sl_vectors[:, :2]
        pmg_substrate_SL = match_pymatgen.substrate_sl_vectors[:, :2]

        # Integer supercell matrices in the original basis (round, then verify by reconstruction)
        M_film_f = np.linalg.solve(film_slab_vectors, pmg_film_SL)
        M_sub_f = np.linalg.solve(substrate_slab_vectors, pmg_substrate_SL)
        film_supercell_matrix = np.rint(M_film_f).astype(int)
        substrate_supercell_matrix = np.rint(M_sub_f).astype(int)

        area = match_pymatgen.match_area
        if self.reduce_result_cell:
            g = get_global_gcd(film_supercell_matrix, substrate_supercell_matrix)
            if g > 1:
                film_supercell_matrix = film_supercell_matrix // g
                substrate_supercell_matrix = substrate_supercell_matrix // g
                area = area / (g**2)

        # Reconstruct matched supercell vectors
        film_SL = film_supercell_matrix @ film_slab_vectors
        sub_SL = substrate_supercell_matrix @ substrate_slab_vectors
        if abs(np.linalg.det(film_SL)) < self.math_precision or abs(np.linalg.det(sub_SL)) < self.math_precision:
            return None

        # Deformation gradient mapping film_SL → sub_SL  (film_SL @ F = sub_SL)
        F = np.linalg.solve(film_SL, sub_SL)

        # For ranking only (rotation-free strain): E = 0.5*(F^T F - I)
        C = F.T @ F
        E2 = 0.5 * (C - np.eye(2))

        # Pack to 3×3 for downstream code
        F3 = convert_2x2_to_3x3(F.tolist())
        E3 = convert_2x2_to_3x3(E2.tolist())
        strain_percentage = self.calculate_total_strain_percentage(E3)

        return ZSLMatchHolder(
            match_id=match_id,
            substrate_supercell_matrix=SupercellMatrix2DSchema(root=substrate_supercell_matrix.tolist()),
            film_supercell_matrix=SupercellMatrix2DSchema(root=film_supercell_matrix.tolist()),
            # IMPORTANT: return the full affine map so the film is rotated + stretched to match
            strain_transformation_matrix=Matrix3x3Schema(root=F3),
            match_area=area,
            total_strain_percentage=strain_percentage,
        )

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
            film_strain_matrix=match_holder.strain_transformation_matrix,
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

    def sort_by_strain_then_area(self, match_holders: List[ZSLMatchHolder], strain_tol=1e-3):
        same_strain_groups = defaultdict(list)
        for m in match_holders:
            strain_key = round(m.total_strain_percentage / strain_tol) * strain_tol
            same_strain_groups[strain_key].append(m)

        result = []
        for strain_key in sorted(same_strain_groups.keys()):
            bucket = same_strain_groups[strain_key]
            bucket.sort(key=lambda m: m.match_area)
            result.extend(bucket)

        return result
