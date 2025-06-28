from typing import List, Tuple

import numpy as np
from mat3ra.code.entity import InMemoryEntityPydantic
from mat3ra.esse.models.core.abstract.matrix_3x3 import Matrix3x3Schema
from mat3ra.esse.models.materials_category_components.entities.auxiliary.two_dimensional.supercell_matrix_2d import (
    SupercellMatrix2DSchema,
)
from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.interface.zsl import ZSLInterfaceAnalyzer
from mat3ra.made.tools.analyze.interface.utils.holders import MatchedSubstrateFilmConfigurationHolder
from mat3ra.made.tools.build.slab.configurations import SlabConfiguration


class GrainBoundaryMatchHolder(InMemoryEntityPydantic):
    """Holder for grain boundary match information."""
    match_id: int
    substrate_transformation_matrix: SupercellMatrix2DSchema
    film_transformation_matrix: SupercellMatrix2DSchema
    match_area: float
    strain_transformation_matrix: Matrix3x3Schema
    total_strain_percentage: float


class GrainBoundaryAnalyzer(InMemoryEntityPydantic):
    """
    Analyzer for creating grain boundaries between two orientations of the same material.

    Uses ZSLInterfaceAnalyzer to find matching supercells between two phases with different
    orientations, then stacks them in x/y direction to create a grain boundary configuration.
    """

    phase_1_material: Material
    phase_2_material: Material
    phase_1_miller_indices: Tuple[int, int, int]
    phase_2_miller_indices: Tuple[int, int, int]
    phase_1_thickness: int = 1
    phase_2_thickness: int = 1
    max_area: float = 50.0
    max_area_ratio_tol: float = 0.09
    max_length_tol: float = 0.03
    max_angle_tol: float = 0.01

    @property
    def phase_1_slab_configuration(self) -> SlabConfiguration:
        """Get slab configuration for phase 1."""
        return SlabConfiguration.from_parameters(
            material_or_dict=self.phase_1_material,
            miller_indices=self.phase_1_miller_indices,
            number_of_layers=self.phase_1_thickness,
            vacuum=0.0,
        )

    @property
    def phase_2_slab_configuration(self) -> SlabConfiguration:
        """Get slab configuration for phase 2."""
        return SlabConfiguration.from_parameters(
            material_or_dict=self.phase_2_material,
            miller_indices=self.phase_2_miller_indices,
            number_of_layers=self.phase_2_thickness,
            vacuum=0.0,
        )

    @property
    def zsl_analyzer(self) -> ZSLInterfaceAnalyzer:
        """Get the underlying ZSL interface analyzer."""
        return ZSLInterfaceAnalyzer(
            substrate_slab_configuration=self.phase_1_slab_configuration,
            film_slab_configuration=self.phase_2_slab_configuration,
            max_area=self.max_area,
            max_area_ratio_tol=self.max_area_ratio_tol,
            max_length_tol=self.max_length_tol,
            max_angle_tol=self.max_angle_tol,
        )

    @property
    def grain_boundary_match_holders(self) -> List[GrainBoundaryMatchHolder]:
        """Get grain boundary match holders from ZSL matches."""
        zsl_matches = self.zsl_analyzer.zsl_match_holders
        match_holders = []

        for idx, zsl_match in enumerate(zsl_matches):
            match_holder = GrainBoundaryMatchHolder(
                match_id=idx,
                substrate_transformation_matrix=zsl_match.substrate_transformation_matrix,
                film_transformation_matrix=zsl_match.film_transformation_matrix,
                match_area=zsl_match.match_area,
                strain_transformation_matrix=zsl_match.strain_transformation_matrix,
                total_strain_percentage=zsl_match.total_strain_percentage,
            )
            match_holders.append(match_holder)

        return match_holders

    def get_grain_boundary_configuration_by_match_id(self, match_id: int) -> MatchedSubstrateFilmConfigurationHolder:
        """Get grain boundary configuration for a specific match ID."""
        match_holders = self.grain_boundary_match_holders
        if match_id < 0 or match_id >= len(match_holders):
            raise ValueError(f"Match ID {match_id} out of range. Available IDs: 0-{len(match_holders)-1}")

        return self.zsl_analyzer.get_strained_configuration_by_match_id(match_id)

    def get_grain_boundary_configurations(self) -> List[MatchedSubstrateFilmConfigurationHolder]:
        """Get all grain boundary configurations."""
        return self.zsl_analyzer.get_strained_configurations()
