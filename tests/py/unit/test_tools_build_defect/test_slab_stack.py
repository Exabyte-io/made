import pytest
from mat3ra.made.tools.build import MaterialWithBuildMetadata
from mat3ra.made.tools.build_components.entities.reusable.two_dimensional.slab_stack.helpers import (
    create_slab_stack,
    recreate_slab_with_fractional_layers,
)
from unit.fixtures.slab import SI_CONVENTIONAL_SLAB_001, SI_SLAB_001_ADDED_FRACTIONAL_LAYER, SI_SLAB_001_ADDED_LAYER
from unit.utils import assert_two_entities_deep_almost_equal


@pytest.mark.parametrize(
    "original_slab_config, number_of_layers, expected_slab_config",
    [
        (SI_CONVENTIONAL_SLAB_001, 1, SI_SLAB_001_ADDED_LAYER),
        (
            SI_CONVENTIONAL_SLAB_001,
            1.5,
            SI_SLAB_001_ADDED_FRACTIONAL_LAYER,
        ),
    ],
)
def test_create_slab_stack(original_slab_config, number_of_layers, expected_slab_config):
    original_slab = MaterialWithBuildMetadata.create(original_slab_config)

    slab_component = recreate_slab_with_fractional_layers(original_slab, number_of_layers)
    # Here we can create isolated defect out of slab_component

    slab_stack = create_slab_stack(original_slab, slab_component)

    assert_two_entities_deep_almost_equal(slab_stack, expected_slab_config, atol=1e-6)
