import pytest
from mat3ra.esse.models.core.abstract.matrix_3x3 import Matrix3x3Schema
from mat3ra.made.material import Material
from mat3ra.made.tools.operations.core.unary import strain
from unit.fixtures.bulk import BULK_Si_CONVENTIONAL
from unit.fixtures.strain import BULK_Si_CONVENTIONAL_STRAINED
from unit.utils import assert_two_entities_deep_almost_equal

STRAIN_TEST_CASES = [
    (BULK_Si_CONVENTIONAL, [[1.1, 0, 0], [0, 1.1, 0], [0, 0, 1.0]], BULK_Si_CONVENTIONAL_STRAINED),
]


@pytest.mark.parametrize("material_config, strain_matrix, expected_material_config", STRAIN_TEST_CASES)
def test_strain(material_config, strain_matrix, expected_material_config):
    material = Material.create(material_config)
    strain_matrix_obj = Matrix3x3Schema(root=strain_matrix)
    expected_material = Material.create(expected_material_config)

    strained_material = strain(material, strain_matrix_obj)

    assert_two_entities_deep_almost_equal(strained_material, expected_material)
