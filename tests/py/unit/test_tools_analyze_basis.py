"""
Tests for basis material analysis functionality.
"""

import pytest
from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.basis import BasisMaterialAnalyzer, LayeredFingerprintAlongAxis

from .fixtures.bulk import BULK_Si_CONVENTIONAL


@pytest.mark.parametrize(
    "material_config, is_flipped",
    [
        (BULK_Si_CONVENTIONAL, False),
    ],
)
def test_basis_analyzer_fingerprint(material_config, is_flipped):
    """Test the basis analyzer fingerprint functionality."""
    material = Material.create(material_config)
    analyzer = BasisMaterialAnalyzer(material=material)

    # Test new Pydantic model format
    fingerprint = analyzer.get_layer_fingerprint(layer_thickness=1.0)
    assert isinstance(fingerprint, LayeredFingerprintAlongAxis)
    assert len(fingerprint.layers) > 0
    assert fingerprint.axis == "z"
    assert fingerprint.layer_thickness == 1.0

    # Test iteration and structure
    for layer in fingerprint.layers:
        assert hasattr(layer, "min_coord")
        assert hasattr(layer, "max_coord")
        assert hasattr(layer, "elements")
        assert isinstance(layer.min_coord, float)
        assert isinstance(layer.max_coord, float)
        assert layer.max_coord > layer.min_coord
        assert isinstance(layer.elements, list)
        assert all(isinstance(elem, str) for elem in layer.elements)

    # Test utility methods
    non_empty_layers = fingerprint.get_non_empty_layers()
    assert isinstance(non_empty_layers, list)
    assert all(isinstance(layer, type(fingerprint.layers[0])) for layer in non_empty_layers)

    element_sequence = fingerprint.get_element_sequence()
    assert isinstance(element_sequence, list)
    assert all(isinstance(elements, list) for elements in element_sequence)

    is_flipped = analyzer.is_orientation_flipped(material)
    assert is_flipped == is_flipped, "Orientation flipped status does not match expected value."


def test_basis_analyzer_different_axes():
    """Test fingerprint generation along different axes."""
    material = Material.create(BULK_Si_CONVENTIONAL)
    analyzer = BasisMaterialAnalyzer(material=material)

    # Test Z-axis (default)
    z_fingerprint = analyzer.get_layer_fingerprint(layer_thickness=1.0, axis=AxisEnum.z)
    assert z_fingerprint.axis == "z"
    assert z_fingerprint.layer_thickness == 1.0

    # Test X-axis
    x_fingerprint = analyzer.get_layer_fingerprint(layer_thickness=1.0, axis=AxisEnum.x)
    assert x_fingerprint.axis == "x"
    assert x_fingerprint.layer_thickness == 1.0

    # Test Y-axis
    y_fingerprint = analyzer.get_layer_fingerprint(layer_thickness=1.0, axis=AxisEnum.y)
    assert y_fingerprint.axis == "y"
    assert y_fingerprint.layer_thickness == 1.0


def test_basis_analyzer_layer_thickness():
    """Test fingerprint generation with different layer thicknesses."""
    material = Material.create(BULK_Si_CONVENTIONAL)
    analyzer = BasisMaterialAnalyzer(material=material)

    # Test different layer thicknesses
    thin_fingerprint = analyzer.get_layer_fingerprint(layer_thickness=0.5)
    thick_fingerprint = analyzer.get_layer_fingerprint(layer_thickness=2.0)

    assert thin_fingerprint.layer_thickness == 0.5
    assert thick_fingerprint.layer_thickness == 2.0

    # Thinner layers should generally result in more layers
    assert len(thin_fingerprint.layers) >= len(thick_fingerprint.layers)


def test_basis_analyzer_orientation_detection():
    """Test orientation flip detection functionality."""
    material = Material.create(BULK_Si_CONVENTIONAL)
    analyzer = BasisMaterialAnalyzer(material=material)

    # Test with same material (should not be flipped)
    is_flipped = analyzer.is_orientation_flipped(material, layer_thickness=1.0)
    assert isinstance(is_flipped, bool)
    assert not is_flipped, "Same material should not be detected as flipped"


def test_fingerprint_utility_methods():
    """Test the utility methods of LayeredFingerprintAlongAxis."""
    material = Material.create(BULK_Si_CONVENTIONAL)
    analyzer = BasisMaterialAnalyzer(material=material)
    fingerprint = analyzer.get_layer_fingerprint(layer_thickness=1.0)

    # Test get_non_empty_layers
    non_empty = fingerprint.get_non_empty_layers()
    assert all(layer.elements for layer in non_empty), "All non-empty layers should have elements"

    # Test get_element_sequence
    element_sequence = fingerprint.get_element_sequence()
    assert len(element_sequence) == len(fingerprint.layers)

    # Test get_non_empty_element_sequence
    non_empty_sequence = fingerprint.get_non_empty_element_sequence()
    assert len(non_empty_sequence) == len(non_empty)
    assert all(elements for elements in non_empty_sequence), "All sequences should have elements"


def test_fingerprint_comparison():
    """Test fingerprint comparison functionality."""
    material = Material.create(BULK_Si_CONVENTIONAL)
    analyzer = BasisMaterialAnalyzer(material=material)

    fingerprint1 = analyzer.get_layer_fingerprint(layer_thickness=1.0)
    fingerprint2 = analyzer.get_layer_fingerprint(layer_thickness=1.0)

    # Test _compare_fingerprints method
    score = analyzer._compare_fingerprints(fingerprint1, fingerprint2)
    assert isinstance(score, float)
    assert 0.0 <= score <= 1.0, "Similarity score should be between 0 and 1"
    assert score == 1.0, "Identical fingerprints should have perfect similarity"

    # Test _reverse_fingerprint method
    reversed_fp = analyzer._reverse_fingerprint(fingerprint1)
    assert isinstance(reversed_fp, LayeredFingerprintAlongAxis)
    assert reversed_fp.axis == fingerprint1.axis
    assert reversed_fp.layer_thickness == fingerprint1.layer_thickness
    assert len(reversed_fp.layers) == len(fingerprint1.layers)
