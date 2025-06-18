from functools import cached_property

import numpy as np
from mat3ra.code.entity import InMemoryEntityPydantic
from mat3ra.esse.models.core.abstract.matrix_3x3 import Matrix3x3Schema
from mat3ra.esse.models.material.reusable.supercell_matrix_2d import SupercellMatrix2DSchema
from mat3ra.made.tools.build.slab.builders import SlabBuilder
from mat3ra.made.tools.build.slab.configuration import SlabConfiguration, SlabStrainedSupercellConfiguration


class InterfaceAnalyzer(InMemoryEntityPydantic):
    substrate_slab_configuration: SlabConfiguration
    film_slab_configuration: SlabConfiguration

    @cached_property
    def substrate_material(self):
        return SlabBuilder().get_material(self.substrate_slab_configuration)

    @cached_property
    def film_material(self):
        return SlabBuilder().get_material(self.film_slab_configuration)

    @cached_property
    def film_strain_matrix(self) -> Matrix3x3Schema:
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

    @property
    def substrate_supercell_matrix(self) -> SupercellMatrix2DSchema:
        return SupercellMatrix2DSchema(root=[[1.0, 0.0], [0.0, 1.0]])

    @property
    def film_supercell_matrix(self) -> SupercellMatrix2DSchema:
        return SupercellMatrix2DSchema(root=[[1.0, 0.0], [0.0, 1.0]])

    @cached_property
    def substrate_strained_configuration(self) -> SlabStrainedSupercellConfiguration:
        return SlabStrainedSupercellConfiguration(
            **self.substrate_slab_configuration.to_dict(),
            xy_supercell_matrix=self.substrate_supercell_matrix,
            strain_matrix=self.substrate_strain_matrix,
        )

    @cached_property
    def film_strained_configuration(self) -> SlabStrainedSupercellConfiguration:
        return SlabStrainedSupercellConfiguration(
            **self.film_slab_configuration.to_dict(),
            xy_supercell_matrix=self.film_supercell_matrix,
            strain_matrix=self.film_strain_matrix,
        )
