import numpy as np

from mat3ra.esse.models.core.abstract.matrix_3x3 import Matrix3x3Schema
from mat3ra.esse.models.material.reusable.supercell_matrix_2d import SupercellMatrix2DSchema
from mat3ra.made.tools.build.interface.configuration import SlabStrainedSupercellConfiguration
from mat3ra.made.tools.build.slab.configuration import SlabConfiguration


class InterfaceAnalyzer:
    """
    Analyzes the required strain to match a film slab to a substrate slab.
    """

    def __init__(self, substrate_config: SlabConfiguration, film_config: SlabConfiguration):
        """
        Initializes the analyzer with substrate and film configurations.

        Args:
            substrate_config (SlabConfiguration): The configuration of the substrate slab.
            film_config (SlabConfiguration): The configuration of the film slab.
        """
        self.substrate_config = substrate_config
        self.film_config = film_config

    def get_strained_configurations(
        self,
    ) -> tuple[SlabStrainedSupercellConfiguration, SlabStrainedSupercellConfiguration]:
        """
        Calculates the strain for the film to match the substrate's lattice and returns
        the configurations for both slabs.

        The substrate is assumed to be rigid, and the film is strained to match the
        substrate's in-plane lattice vectors. A 1x1 supercell is assumed for both.

        Returns:
            tuple[SlabStrainedSupercellConfiguration, SlabStrainedSupercellConfiguration]:
                A tuple containing the strained supercell configurations for the
                substrate and the film.
        """
        substrate_material = self.substrate_config.atomic_layers.crystal
        film_material = self.film_config.atomic_layers.crystal

        substrate_vectors = np.array(substrate_material.lattice.vector_arrays)
        film_vectors = np.array(film_material.lattice.vector_arrays)

        # Assuming in-plane vectors are the first two lattice vectors
        substrate_2d_vectors = substrate_vectors[:2, :2]
        film_2d_vectors = film_vectors[:2, :2]

        try:
            inv_film_2d_vectors = np.linalg.inv(film_2d_vectors)
        except np.linalg.LinAlgError:
            raise ValueError("Film lattice vectors are not linearly independent.")

        strain_2d_matrix = inv_film_2d_vectors @ substrate_2d_vectors

        film_strain_matrix_3x3 = np.identity(3)
        film_strain_matrix_3x3[:2, :2] = strain_2d_matrix

        substrate_strain_matrix_3x3 = np.identity(3)

        identity_supercell_matrix = [[1.0, 0.0], [0.0, 1.0]]

        substrate_strained_config = SlabStrainedSupercellConfiguration(
            **self.substrate_config.model_dump(),
            xy_supercell_matrix=SupercellMatrix2DSchema(root=identity_supercell_matrix),
            strain_matrix=Matrix3x3Schema(root=substrate_strain_matrix_3x3.tolist()),
        )

        film_strained_config = SlabStrainedSupercellConfiguration(
            **self.film_config.model_dump(),
            xy_supercell_matrix=SupercellMatrix2DSchema(root=identity_supercell_matrix),
            strain_matrix=Matrix3x3Schema(root=film_strain_matrix_3x3.tolist()),
        )

        return substrate_strained_config, film_strained_config
