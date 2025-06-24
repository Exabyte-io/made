from functools import cached_property
from typing import List, Optional

import numpy as np
from mat3ra.code.entity import InMemoryEntityPydantic
from mat3ra.esse.models.materials_category_components.entities.auxiliary.two_dimensional.supercell_matrix_2d import (
    SupercellMatrix2DSchema,
)
from mat3ra.made.lattice import Lattice
from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.interface.enums import angle_to_supercell_matrix_values_for_hex
from mat3ra.made.tools.analyze.interface.simple import InterfaceAnalyzer
from mat3ra.made.tools.analyze.interface.utils.holders import MatchedSubstrateFilmConfigurationHolder
from mat3ra.made.utils import create_2d_supercell_matrices, get_angle_from_rotation_matrix_2d


class CommensurateLatticeMatchHolder(InMemoryEntityPydantic):
    match_id: int
    angle: float
    xy_supercell_matrix_film: List[List[int]]
    xy_supercell_matrix_substrate: List[List[int]]


class CommensurateInterfaceAnalyzer(InterfaceAnalyzer):
    """
    Interface analyzer using commensurate lattice matching for twisted interfaces.

    For twisted bilayer materials where film and substrate are the same material,
    this analyzer finds commensurate supercells at a target twist angle.
    """

    target_angle: float = 0.0
    angle_tolerance: float = 0.1
    max_supercell_matrix_int: Optional[int] = None
    limit_max_int: int = 50
    return_first_match: bool = True

    def __init__(self, **data):
        # For twisted bilayers, both film and substrate use the same material
        if "film_slab_configuration" not in data and "substrate_slab_configuration" in data:
            data["film_slab_configuration"] = data["substrate_slab_configuration"]
        elif "substrate_slab_configuration" not in data and "film_slab_configuration" in data:
            data["substrate_slab_configuration"] = data["film_slab_configuration"]
        super().__init__(**data)

    @property
    def material(self) -> Material:
        """For twisted bilayers, both layers use the same base material."""
        return self.substrate_material

    def _get_initial_guess_for_max_int(self, material: Material, target_angle: float) -> int:
        """
        Determine the maximum integer for the transformation matrices based on the target angle.
        """
        if material.lattice.type == Lattice.__types__.HEX:
            xy_supercell_matrix_for_closest_angle = min(
                angle_to_supercell_matrix_values_for_hex, key=lambda x: abs(x["angle"] - target_angle)
            )
            return max(abs(x) for row in xy_supercell_matrix_for_closest_angle["xy_supercell"] for x in row)
        return 1

    def _generate_commensurate_lattices(
        self,
        material: Material,
        target_angle: float,
        max_supercell_matrix_element_int: int,
    ) -> List[CommensurateLatticeMatchHolder]:
        """
        Generate all commensurate lattices for a given target angle.
        """
        a = material.lattice.vector_arrays[0][:2]
        b = material.lattice.vector_arrays[1][:2]

        matrices = create_2d_supercell_matrices(max_supercell_matrix_element_int)
        matrix_ab = np.array([a, b])
        matrix_ab_inverse = np.linalg.inv(matrix_ab)

        solutions: List[CommensurateLatticeMatchHolder] = []
        match_id = 0

        for index1, matrix1 in enumerate(matrices):
            for index2, matrix2 in enumerate(matrices[0 : index1 + 1]):
                matrix2_inverse = np.linalg.inv(matrix2)
                intermediate_product = matrix2_inverse @ matrix1
                product = matrix_ab_inverse @ intermediate_product @ matrix_ab
                angle = get_angle_from_rotation_matrix_2d(product)

                if angle is not None and np.abs(angle - target_angle) < self.angle_tolerance:
                    solutions.append(
                        CommensurateLatticeMatchHolder(
                            match_id=match_id,
                            angle=angle,
                            xy_supercell_matrix_substrate=matrix1.tolist(),
                            xy_supercell_matrix_film=matrix2.tolist(),
                        )
                    )
                    match_id += 1

                    if self.return_first_match:
                        return solutions

        return solutions

    @cached_property
    def commensurate_match_holders(self) -> List[CommensurateLatticeMatchHolder]:
        """Find all commensurate lattice matches for the target angle."""
        max_int = self.max_supercell_matrix_int or self._get_initial_guess_for_max_int(self.material, self.target_angle)

        commensurate_lattice_matches: List[CommensurateLatticeMatchHolder] = []

        while not commensurate_lattice_matches and max_int < self.limit_max_int:
            commensurate_lattice_matches = self._generate_commensurate_lattices(
                self.material, self.target_angle, max_int
            )
            max_int += 1

        return commensurate_lattice_matches

    def get_strained_configuration_by_match_id(self, match_id: int) -> MatchedSubstrateFilmConfigurationHolder:
        """Get strained slab configurations for a specific match ID."""
        match_holders = self.commensurate_match_holders
        if match_id < 0 or match_id >= len(match_holders):
            raise ValueError(f"Match ID {match_id} out of range. Available IDs: 0-{len(match_holders)-1}")

        match_holder = match_holders[match_id]
        return self._create_strained_configs_from_match(match_holder)

    def _create_strained_configs_from_match(
        self, match_holder: CommensurateLatticeMatchHolder
    ) -> MatchedSubstrateFilmConfigurationHolder:
        substrate_matrix_schema = SupercellMatrix2DSchema(root=match_holder.xy_supercell_matrix_substrate)
        film_matrix_schema = SupercellMatrix2DSchema(root=match_holder.xy_supercell_matrix_film)

        substrate_config = self.get_component_strained_configuration(
            self.substrate_slab_configuration, self._no_strain_matrix, xy_supercell_matrix=substrate_matrix_schema
        )

        film_config = self.get_component_strained_configuration(
            self.film_slab_configuration, self._no_strain_matrix, xy_supercell_matrix=film_matrix_schema
        )

        return MatchedSubstrateFilmConfigurationHolder(
            match_id=match_holder.match_id,
            substrate_configuration=substrate_config,
            film_configuration=film_config,
        )

    def get_strained_configurations(self) -> List[MatchedSubstrateFilmConfigurationHolder]:
        """Get all strained configurations for all matches."""
        strained_configs = []

        for match_holder in self.commensurate_match_holders:
            config_holder = self._create_strained_configs_from_match(match_holder)
            strained_configs.append(config_holder)

        return strained_configs
