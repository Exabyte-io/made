from functools import cached_property
from typing import List

import numpy as np
from mat3ra.code.entity import InMemoryEntityPydantic
from mat3ra.esse.models.core.abstract.matrix_3x3 import Matrix3x3Schema
from mat3ra.esse.models.materials_category_components.entities.auxiliary.two_dimensional.supercell_matrix_2d import (
    SupercellMatrix2DSchema,
)
from mat3ra.utils.matrix import convert_2x2_to_3x3
from pymatgen.analysis.interfaces.coherent_interfaces import CoherentInterfaceBuilder, ZSLGenerator

from mat3ra.made.tools.analyze.interface.simple import InterfaceAnalyzer
from mat3ra.made.tools.analyze.interface.utils.holders import MatchedSubstrateFilmConfigurationHolder
from mat3ra.made.tools.analyze.lattice import get_conventional_material, get_primitive_material
from mat3ra.made.tools.build.slab.builders import SlabBuilder, SlabStrainedSupercellBuilder
from mat3ra.made.tools.convert import to_pymatgen
from mat3ra.made.tools.operations.core.unary import supercell
from mat3ra.made.utils import calculate_von_mises_strain


class ZSLMatchHolder(InMemoryEntityPydantic):
    match_id: int
    substrate_transformation_matrix: SupercellMatrix2DSchema
    film_transformation_matrix: SupercellMatrix2DSchema
    strain_transformation_matrix: Matrix3x3Schema
    substrate_sl_vector_arrays: np.ndarray = None
    film_sl_vector_arrays: np.ndarray = None
    match_area: float
    total_strain_percentage: float


class ZSLInterfaceAnalyzer(InterfaceAnalyzer):
    """Interface analyzer using Pymatgen's ZSL algorithm to find matching supercells."""

    max_area: float = 50.0
    max_area_ratio_tol: float = 0.09
    max_length_tol: float = 0.03
    max_angle_tol: float = 0.01

    @cached_property
    def _pymatgen_zsl_generator(self) -> ZSLGenerator:
        return ZSLGenerator(
            max_area=self.max_area,
            max_area_ratio_tol=self.max_area_ratio_tol,
            max_length_tol=self.max_length_tol,
            max_angle_tol=self.max_angle_tol,
        )

    @cached_property
    def _pymatgen_coherent_interface_builder(self) -> CoherentInterfaceBuilder:
        return CoherentInterfaceBuilder(
            substrate_structure=to_pymatgen(self.substrate_slab_configuration.atomic_layers.crystal),
            film_structure=to_pymatgen(self.film_slab_configuration.atomic_layers.crystal),
            substrate_miller=self.substrate_slab_configuration.atomic_layers.miller_indices,
            film_miller=self.film_slab_configuration.atomic_layers.miller_indices,
            zslgen=self._pymatgen_zsl_generator,
        )

    @classmethod
    def calculate_total_strain_percentage(cls, strain_matrix: np.ndarray) -> float:
        return calculate_von_mises_strain(strain_matrix)

    @cached_property
    def zsl_match_holders(self) -> List[ZSLMatchHolder]:
        zsl_matches_pymatgen = self._pymatgen_coherent_interface_builder.zsl_matches
        match_holders = []

        for idx, match_pymatgen in enumerate(zsl_matches_pymatgen):
            match_holder = ZSLMatchHolder(
                match_id=idx,
                substrate_transformation_matrix=SupercellMatrix2DSchema(
                    root=match_pymatgen.substrate_transformation.tolist()
                ),
                film_transformation_matrix=SupercellMatrix2DSchema(root=match_pymatgen.film_transformation.tolist()),
                strain_transformation_matrix=Matrix3x3Schema(root=match_pymatgen.match_transformation.tolist()),
                substrate_sl_vector_arrays=match_pymatgen.substrate_sl_vectors,
                film_sl_vector_arrays=match_pymatgen.film_sl_vectors,
                match_area=match_pymatgen.match_area,
                total_strain_percentage=self.calculate_total_strain_percentage(match_pymatgen.match_transformation),
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

        film_slab = SlabBuilder().get_material(self.film_slab_configuration)
        substrate_slab = SlabBuilder().get_material(self.substrate_slab_configuration)

        # we need to calculate supercell matrices to transform from our slabs to SL slabs from ZSL matches
        film_current_lattice_vectors = film_slab.lattice.vector_arrays
        substrate_current_lattice_vectors = substrate_slab.lattice.vector_arrays

        film_target_lattice_vectors = match_holder.film_sl_vector_arrays
        substrate_target_lattice_vectors = match_holder.substrate_sl_vector_arrays

        film_float_matrix = (
            np.linalg.inv(np.array(film_current_lattice_vectors)[:2, :2])
            @ np.array(film_target_lattice_vectors)[:2, :2]
        )
        film_transformation_matrix = np.rint(film_float_matrix).astype(int)

        substrate_float_matrix = (
            np.linalg.inv(np.array(substrate_current_lattice_vectors)[:2, :2])
            @ np.array(substrate_target_lattice_vectors)[:2, :2]
        )
        substrate_transformation_matrix = np.rint(substrate_float_matrix).astype(int)

        # 1) build the supercells
        film_sc = supercell(film_slab, convert_2x2_to_3x3(film_transformation_matrix))
        substrate_sc = supercell(substrate_slab, convert_2x2_to_3x3(substrate_transformation_matrix))

        # 2) pull off their 2×2 in-plane lattice matrices (row-vectors)
        C_film_sc = np.array(film_sc.lattice.vector_arrays)[:2, :2]
        C_substrate_sc = np.array(substrate_sc.lattice.vector_arrays)[:2, :2]

        # 3) compute the exact 2×2 deformation gradient
        F2 = np.linalg.inv(C_film_sc) @ C_substrate_sc

        # 4) embed into a full 3×3 strain
        strain_3x3 = np.eye(3)
        strain_3x3[:2, :2] = F2

        film_strain_matrix = Matrix3x3Schema(root=strain_3x3.tolist())

        return super().create_matched_configuration_holder(
            substrate_slab_config=self.substrate_slab_configuration,
            film_slab_config=self.film_slab_configuration,
            match_id=match_holder.match_id,
            substrate_xy_supercell_matrix=SupercellMatrix2DSchema(root=substrate_transformation_matrix.tolist()),
            film_xy_supercell_matrix=SupercellMatrix2DSchema(root=film_transformation_matrix.tolist()),
            substrate_strain_matrix=self._no_strain_matrix,
            film_strain_matrix=Matrix3x3Schema(root=film_strain_matrix.root),
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
