import pytest
from mat3ra.esse.models.core.abstract.matrix_3x3 import Matrix3x3Schema
from mat3ra.made.material import Material
from mat3ra.made.tools.build_components.operations.core.modifications.strain.helpers import create_strain
from unit.fixtures.strain import BULK_Si_CONVENTIONAL_STRAINED
from unit.utils import assert_two_entities_deep_almost_equal

from .fixtures.bulk import BULK_Si_CONVENTIONAL


@pytest.mark.parametrize(
    "material_config, strain_matrix, expected_material_config",
    [
        (BULK_Si_CONVENTIONAL, [[1.1, 0, 0], [0, 1.1, 0], [0, 0, 1.0]], BULK_Si_CONVENTIONAL_STRAINED),
    ],
)
def test_create_strain(material_config, strain_matrix, expected_material_config):
    material = Material.create(material_config)
    strained_material = create_strain(material, Matrix3x3Schema(root=strain_matrix))

    assert_two_entities_deep_almost_equal(strained_material, expected_material_config)

    build_step = strained_material.metadata.build[-1]
    assert build_step.configuration["type"] == "StrainConfiguration"
    assert build_step.configuration["strain_matrix"] == strain_matrix
    assert build_step.build_parameters == {}


def test_create_strain_with_scale_factor():
    material = Material.create(BULK_Si_CONVENTIONAL)
    strained_material = create_strain(material, scale_factor=1.1)

    assert strained_material.name == f"{material.name} (scale=1.1000)"
    build_step = strained_material.metadata.build[-1]
    assert build_step.configuration["scale_factor"] == 1.1
    assert build_step.configuration["strain_matrix"] == [[1.1, 0.0, 0.0], [0.0, 1.1, 0.0], [0.0, 0.0, 1.1]]
