import pytest
from mat3ra.made.material import Material
from mat3ra.made.tools.build.defect.builders import SlabDefectBuilder, SlabDefectBuilderParameters
from mat3ra.utils import assertion as assertion_utils
from unit.fixtures.slab import SI_CONVENTIONAL_SLAB_001, SI_SLAB_001_ADDED_FRACTIONAL_LAYER, SI_SLAB_001_ADDED_LAYER


@pytest.mark.skip(reason="Slab with additional layers should be adjusted")
@pytest.mark.parametrize(
    "original_slab_config, layers_to_add, builder_params_dict, expected_slab_config",
    [(SI_CONVENTIONAL_SLAB_001, 1, {"auto_add_vacuum": True, "vacuum_thickness": 5.0}, SI_SLAB_001_ADDED_LAYER)],
)
def test_create_material_with_additional_layers(
    original_slab_config, layers_to_add, builder_params_dict, expected_slab_config
):
    """Test adding layers to a slab material"""
    # Create the builder
    builder_params = SlabDefectBuilderParameters(**builder_params_dict)
    builder = SlabDefectBuilder(build_parameters=builder_params)

    # Test adding 1 layer to SI_SLAB_001
    original_slab = Material.create(original_slab_config)
    expected_slab = Material.create(expected_slab_config)
    slab_with_additional_layer = builder.create_material_with_additional_layers(original_slab, layers_to_add)

    assertion_utils.assert_deep_almost_equal(slab_with_additional_layer, expected_slab)


@pytest.mark.skip(reason="Slab with additional layers should be adjusted")
@pytest.mark.parametrize(
    "original_slab_config, layers_to_add, builder_params_dict, expected_slab_config",
    [
        (
            SI_CONVENTIONAL_SLAB_001,
            1.5,
            {"auto_add_vacuum": True, "vacuum_thickness": 5.0},
            SI_SLAB_001_ADDED_FRACTIONAL_LAYER,
        )
    ],
)
def test_create_material_with_additional_fractional_layers(
    original_slab_config, layers_to_add, builder_params_dict, expected_slab_config
):
    """Test adding fractional layers to a slab material"""
    # Create the builder
    builder_params = SlabDefectBuilderParameters(**builder_params_dict)
    builder = SlabDefectBuilder(build_parameters=builder_params)

    # Test adding 1.5 layers to SI_SLAB_001
    original_slab = Material.create(original_slab_config)
    expected_slab = Material.create(expected_slab_config)
    slab_with_fractional_layer = builder.create_material_with_additional_layers(original_slab, layers_to_add)

    # Compare with expected fixture
    assertion_utils.assert_deep_almost_equal(slab_with_fractional_layer, expected_slab)
