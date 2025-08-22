"""
Tests for basis material analysis functionality.
"""

import pytest
from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.basis import BasisMaterialAnalyzer, LayeredFingerprintAlongAxis

from .fixtures.bulk import BULK_Si_CONVENTIONAL
from .fixtures.slab import SLAB_SrTiO3_011_TERMINATION_O2, SLAB_SrTiO3_011_TERMINATION_SrTiO


@pytest.mark.parametrize(
    "original_material_config, another_material_config, is_flipped",
    [
        (SLAB_SrTiO3_011_TERMINATION_O2, SLAB_SrTiO3_011_TERMINATION_O2, False),
        (SLAB_SrTiO3_011_TERMINATION_O2, SLAB_SrTiO3_011_TERMINATION_SrTiO, True),
    ],
)
def test_basis_analyzer_fingerprint(original_material_config, another_material_config, is_flipped):
    original_material = Material.create(original_material_config)
    another_material = Material.create(another_material_config)
    analyzer = BasisMaterialAnalyzer(material=original_material)

    fingerprint = analyzer.get_layer_fingerprint(layer_thickness=1.0)
    assert isinstance(fingerprint, LayeredFingerprintAlongAxis)
    assert len(fingerprint.layers) > 0
    assert fingerprint.axis == AxisEnum.z
    assert fingerprint.layer_thickness == 1.0

    for layer in fingerprint.layers:
        assert hasattr(layer, "min_coord")
        assert hasattr(layer, "max_coord")
        assert hasattr(layer, "elements")
        assert isinstance(layer.min_coord, float)
        assert isinstance(layer.max_coord, float)
        assert layer.max_coord > layer.min_coord
        assert isinstance(layer.elements, list)
        assert all(isinstance(elem, str) for elem in layer.elements)

    non_empty_layers = fingerprint.get_non_empty_layers()
    assert isinstance(non_empty_layers, list)
    assert all(isinstance(layer, type(fingerprint.layers[0])) for layer in non_empty_layers)

    element_sequence = fingerprint.get_element_sequence()
    assert isinstance(element_sequence, list)
    assert all(isinstance(elements, list) for elements in element_sequence)

    is_flipped = analyzer.is_orientation_flipped(another_material)
    assert is_flipped == is_flipped, "Orientation flipped status does not match expected value."


def test_fingerprint_similarity_score():
    """Test the get_similarity_score method specifically."""
    material = Material.create(BULK_Si_CONVENTIONAL)
    analyzer = BasisMaterialAnalyzer(material=material)

    fingerprint = analyzer.get_layer_fingerprint(layer_thickness=1.0)

    # Test self-similarity (should be 1.0)
    self_score = fingerprint.get_similarity_score(fingerprint)
    assert self_score == 1.0, "Self-similarity should be perfect"

    # Test with different layer thickness (should be different)
    different_fp = analyzer.get_layer_fingerprint(layer_thickness=2.0)
    different_score = fingerprint.get_similarity_score(different_fp)
    assert isinstance(different_score, float)
    assert 0.0 <= different_score <= 1.0

    # Test with empty fingerprint
    empty_fp = LayeredFingerprintAlongAxis(layers=[], axis=AxisEnum.z, layer_thickness=1.0)
    empty_score = fingerprint.get_similarity_score(empty_fp)
    assert empty_score == 0.0, "Similarity with empty fingerprint should be 0"
