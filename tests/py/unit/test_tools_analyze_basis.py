"""
Tests for BasisMaterialAnalyzer class.
"""

from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.basis import BasisMaterialAnalyzer
from mat3ra.made.tools.analyze.fingerprint import LayeredFingerprintAlongAxis

from .fixtures.bulk import BULK_Si_CONVENTIONAL
from .fixtures.slab import SLAB_SrTiO3_011_TERMINATION_O2


def test_basis_analyzer_fingerprint():
    """Test fingerprint creation functionality."""
    material = Material.create(SLAB_SrTiO3_011_TERMINATION_O2)
    analyzer = BasisMaterialAnalyzer(material=material)

    # Test basic fingerprint functionality
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


def test_get_layer_fingerprint():
    """Test BasisMaterialAnalyzer.get_layer_fingerprint method."""
    material = Material.create(BULK_Si_CONVENTIONAL)
    analyzer = BasisMaterialAnalyzer(material=material)

    # Test default z-axis fingerprint
    fingerprint = analyzer.get_layer_fingerprint(layer_thickness=1.0)
    assert isinstance(fingerprint, LayeredFingerprintAlongAxis)
    assert len(fingerprint.layers) > 0
    assert fingerprint.axis == AxisEnum.z
    assert fingerprint.layer_thickness == 1.0

    # Test different axes
    for axis in [AxisEnum.x, AxisEnum.y, AxisEnum.z]:
        fingerprint = analyzer.get_layer_fingerprint(layer_thickness=1.0, axis=axis)
        assert isinstance(fingerprint, LayeredFingerprintAlongAxis)
        assert fingerprint.axis == axis
        assert fingerprint.layer_thickness == 1.0
        assert len(fingerprint.layers) > 0


def test_get_material_fingerprint():
    """Test BasisMaterialAnalyzer.get_material_fingerprint method for all axes."""
    material = Material.create(BULK_Si_CONVENTIONAL)
    analyzer = BasisMaterialAnalyzer(material=material)

    # Test full material fingerprint across all axes
    material_fingerprint = analyzer.get_material_fingerprint(layer_thickness=1.0)
    
    assert hasattr(material_fingerprint, "x_axis")
    assert hasattr(material_fingerprint, "y_axis")
    assert hasattr(material_fingerprint, "z_axis")
    assert material_fingerprint.layer_thickness == 1.0
    
    # Verify each axis fingerprint
    for axis_name in ["x_axis", "y_axis", "z_axis"]:
        axis_fingerprint = getattr(material_fingerprint, axis_name)
        assert isinstance(axis_fingerprint, LayeredFingerprintAlongAxis)
        assert len(axis_fingerprint.layers) > 0
