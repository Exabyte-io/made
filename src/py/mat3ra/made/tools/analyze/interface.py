from functools import cached_property

import numpy as np
from mat3ra.code.entity import InMemoryEntityPydantic
from mat3ra.esse.models.core.abstract.matrix_3x3 import Matrix3x3Schema

from mat3ra.made.tools.build.slab.builders import SlabBuilder
from mat3ra.made.tools.build.slab.configuration import SlabConfiguration, SlabStrainedSupercellConfiguration


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

    @cached_property
    def film_strain_matrix(self) -> Matrix3x3Schema:
        """
        Calculate the strain matrix for the film material to match the substrate material's lattice vectors.
        """
        substrate_vectors = np.array(self.substrate_material.lattice.vector_arrays)
        film_vectors = np.array(self.film_material.lattice.vector_arrays)

        substrate_2d_vectors = substrate_vectors[:2, :2]
        film_2d_vectors = film_vectors[:2, :2]

        try:
            inv_film_2d_vectors = np.linalg.inv(film_2d_vectors)
        except np.linalg.LinAlgError:
            raise ValueError("Film lattice vectors are not linearly independent.")

        strain_2d_matrix = inv_film_2d_vectors @ substrate_2d_vectors

        film_strain_matrix_3x3 = np.identity(3)
        film_strain_matrix_3x3[:2, :2] = strain_2d_matrix
        return Matrix3x3Schema(root=film_strain_matrix_3x3.tolist())

    @property
    def substrate_strain_matrix(self) -> Matrix3x3Schema:
        return Matrix3x3Schema(root=np.identity(3).tolist())

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

    @cached_property
    def substrate_strained_configuration(self) -> SlabStrainedSupercellConfiguration:
        return self.get_component_strained_configuration(
            self.substrate_slab_configuration, self.substrate_strain_matrix
        )

    @cached_property
    def film_strained_configuration(self) -> SlabStrainedSupercellConfiguration:
        return self.get_component_strained_configuration(self.film_slab_configuration, self.film_strain_matrix)
