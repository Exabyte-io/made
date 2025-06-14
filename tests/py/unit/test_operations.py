import numpy as np
from mat3ra.esse.models.core.abstract.matrix_3x3 import Matrix3x3Schema

from mat3ra.made.material import Material
from mat3ra.made.tools.operations.core.unary import strain
from mat3ra.utils import assertion as assertion_utils
from unit.fixtures.bulk import BULK_Si_CONVENTIONAL


def test_strain():
    """
    Tests that the strain operation correctly modifies the lattice vectors of a material.
    """
    material = Material.create(BULK_Si_CONVENTIONAL)
    strain_matrix_data = [[1.01, 0, 0], [0, 1, 0], [0, 0, 1]]
    strain_matrix = Matrix3x3Schema(root=strain_matrix_data)

    strained_material = strain(material, strain_matrix)

    original_lattice_vectors = np.array(material.lattice.vector_arrays)
    expected_strained_lattice_vectors = np.dot(original_lattice_vectors, np.array(strain_matrix_data))

    assertion_utils.assert_deep_almost_equal(
        strained_material.lattice.vector_arrays, expected_strained_lattice_vectors.tolist()
    )
