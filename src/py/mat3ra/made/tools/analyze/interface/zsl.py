from functools import cached_property
from typing import List

import numpy as np
from mat3ra.code.entity import InMemoryEntityPydantic
from mat3ra.esse.models.core.abstract.matrix_3x3 import Matrix3x3Schema
from mat3ra.esse.models.materials_category_components.entities.auxiliary.two_dimensional.supercell_matrix_2d import (
    SupercellMatrix2DSchema,
)
from mat3ra.made.tools.analyze.interface import InterfaceAnalyzer
from mat3ra.made.tools.build.slab.configuration import SlabStrainedSupercellConfiguration
from mat3ra.made.tools.convert import to_pymatgen
from mat3ra.made.tools.operations.core.unary import supercell
from mat3ra.made.tools.utils import supercell_matrix_2d_schema_to_list
from pymatgen.analysis.interfaces.coherent_interfaces import CoherentInterfaceBuilder, ZSLGenerator


class MatchedSubstrateFilmConfigurationHolder(InMemoryEntityPydantic):
    match_id: int
    substrate_configuration: SlabStrainedSupercellConfiguration
    film_configuration: SlabStrainedSupercellConfiguration


class ZSLMatchHolder(InMemoryEntityPydantic):
    match_id: int
    substrate_transformation_matrix: SupercellMatrix2DSchema
    film_transformation_matrix: SupercellMatrix2DSchema
    match_area: float
    strain_transformation_matrix: Matrix3x3Schema
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
        """Calculate von Mises strain from a 2D strain transformation matrix."""
        exx = strain_matrix[0, 0] - 1.0
        eyy = strain_matrix[1, 1] - 1.0
        exy = strain_matrix[0, 1]

        von_mises = np.sqrt(exx**2 - exx * eyy + eyy**2 + 3 * exy**2) / np.sqrt(2)
        return abs(von_mises) * 100.0

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
                match_area=match_pymatgen.match_area,
                strain_transformation_matrix=Matrix3x3Schema(root=match_pymatgen.match_transformation.tolist()),
                total_strain_percentage=self.calculate_total_strain_percentage(match_pymatgen.match_transformation),
            )
            match_holders.append(match_holder)

        return match_holders

    def get_strained_configuration_by_match_id(self, match_id: int) -> MatchedSubstrateFilmConfigurationHolder:
        match_holders = self.zsl_match_holders
        if match_id < 0 or match_id >= len(match_holders):
            raise ValueError(f"Match ID {match_id} out of range. Available IDs: 0-{len(match_holders)-1}")

        match_holder = match_holders[match_id]
        return self._create_strained_configs_from_match(match_holder)

    def _create_strained_configs_from_match(
        self, match_holder: ZSLMatchHolder
    ) -> MatchedSubstrateFilmConfigurationHolder:
        substrate_supercell = supercell(self.substrate_material, match_holder.substrate_transformation_matrix)
        film_supercell = supercell(self.film_material, match_holder.film_transformation_matrix)

        film_strain_3d = self.get_film_strain_matrix(
            substrate_supercell.lattice.vector_arrays, film_supercell.lattice.vector_arrays
        )

        substrate_matrix = supercell_matrix_2d_schema_to_list(match_holder.substrate_transformation_matrix)
        film_matrix = supercell_matrix_2d_schema_to_list(match_holder.film_transformation_matrix)

        substrate_config = SlabStrainedSupercellConfiguration(
            stack_components=self.substrate_slab_configuration.stack_components,
            direction=self.substrate_slab_configuration.direction,
            xy_supercell_matrix=substrate_matrix,
            strain_matrix=self._no_strain_matrix,
        )

        film_config = SlabStrainedSupercellConfiguration(
            stack_components=self.film_slab_configuration.stack_components,
            direction=self.film_slab_configuration.direction,
            xy_supercell_matrix=film_matrix,
            strain_matrix=film_strain_3d,
        )

        return MatchedSubstrateFilmConfigurationHolder(
            match_id=match_holder.match_id,
            substrate_configuration=substrate_config,
            film_configuration=film_config,
        )

    def get_strained_configurations(
        self,
    ) -> List[MatchedSubstrateFilmConfigurationHolder]:
        strained_configs = []

        for match_holder in self.zsl_match_holders:
            config_holder = self._create_strained_configs_from_match(match_holder)
            strained_configs.append(config_holder)

        return strained_configs
