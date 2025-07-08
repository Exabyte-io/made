import pytest

from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.slab import SlabMaterialAnalyzer
from mat3ra.made.tools.build.slab.builders import SlabBuilder, SlabWithAdditionalLayersBuilder
from mat3ra.made.tools.modify import translate_to_z_level
from unit.fixtures.slab import SI_CONVENTIONAL_SLAB_001, SI_SLAB_001_ADDED_FRACTIONAL_LAYER, SI_SLAB_001_ADDED_LAYER
from unit.utils import assert_two_entities_deep_almost_equal


@pytest.mark.parametrize(
    "original_slab_config, layers_to_add, analyzer_params_dict, expected_slab_config",
    [(SI_CONVENTIONAL_SLAB_001, 1, {"vacuum_thickness": 5.0}, SI_SLAB_001_ADDED_LAYER)],
)
def test_analyzer_get_slab_configurations(
    original_slab_config, layers_to_add, analyzer_params_dict, expected_slab_config
):
    original_slab = Material.create(original_slab_config)
    analyzer = SlabMaterialAnalyzer(material=original_slab)
    expected_slab = translate_to_z_level(Material.create(expected_slab_config), "bottom")

    slab_with_additional_layers_config, slab_with_original_layers_config = (
        analyzer.get_slab_with_additional_layers_configurations(
            additional_layers=layers_to_add, vacuum_thickness=analyzer_params_dict["vacuum_thickness"]
        )
    )

    builder_additional = SlabWithAdditionalLayersBuilder()
    builder_original = SlabBuilder()
    slab_with_additional_layers = builder_additional.get_material(slab_with_additional_layers_config)
    slab_with_original_layers_adjusted = builder_original.get_material(slab_with_original_layers_config)

    assert_two_entities_deep_almost_equal(slab_with_additional_layers, expected_slab, atol=1e-6)
    assert_two_entities_deep_almost_equal(
        slab_with_original_layers_adjusted.lattice.vector_arrays, expected_slab.lattice.vector_arrays
    )


@pytest.mark.parametrize(
    "original_slab_config, layers_to_add, analyzer_params_dict, expected_slab_config",
    [
        (
            SI_CONVENTIONAL_SLAB_001,
            1.5,
            {"vacuum_thickness": 5.0},
            SI_SLAB_001_ADDED_FRACTIONAL_LAYER,
        )
    ],
)
def test_analyzer_get_slab_configurations_fractional(
    original_slab_config, layers_to_add, analyzer_params_dict, expected_slab_config
):
    original_slab = Material.create(original_slab_config)
    analyzer = SlabMaterialAnalyzer(material=original_slab)

    slab_with_additional_layers_config, slab_with_original_layers_config = (
        analyzer.get_slab_with_additional_layers_configurations(
            additional_layers=layers_to_add, vacuum_thickness=analyzer_params_dict["vacuum_thickness"]
        )
    )

    slab_with_additional_layers = SlabWithAdditionalLayersBuilder().get_material(slab_with_additional_layers_config)
    slab_with_original_layers_adjusted = SlabBuilder().get_material(slab_with_original_layers_config)
    expected_slab = Material.create(expected_slab_config)

    assert_two_entities_deep_almost_equal(slab_with_additional_layers, expected_slab_config, atol=1e-6)
    assert_two_entities_deep_almost_equal(
        slab_with_original_layers_adjusted.lattice.vector_arrays, expected_slab.lattice.vector_arrays
    )
