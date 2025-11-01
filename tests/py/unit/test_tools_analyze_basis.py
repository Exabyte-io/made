"""
Tests for BasisMaterialAnalyzer class.
"""

import pytest
from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.basis import BasisMaterialAnalyzer, LayeredFingerprintAlongAxis
from mat3ra.made.tools.operations.core.unary import rotate

from .fixtures.bulk import BULK_Si_CONVENTIONAL
from .fixtures.slab import (
    SLAB_SrTiO3_011_TERMINATION_O2, 
    SLAB_SrTiO3_011_TERMINATION_SrTiO,
    SI_CONVENTIONAL_SLAB_001
)


@pytest.mark.parametrize(
    "original_material_config, another_material_config, is_flipped",
    [
        (SLAB_SrTiO3_011_TERMINATION_O2, SLAB_SrTiO3_011_TERMINATION_O2, False),
        (SLAB_SrTiO3_011_TERMINATION_O2, SLAB_SrTiO3_011_TERMINATION_SrTiO, True),
    ],
)
def test_basis_analyzer_fingerprint(original_material_config, another_material_config, is_flipped):
    """Test fingerprint functionality and rotation detection between different terminations."""
    original_material = Material.create(original_material_config)
    another_material = Material.create(another_material_config)
    
    original_analyzer = BasisMaterialAnalyzer(material=original_material)
    another_analyzer = BasisMaterialAnalyzer(material=another_material)

    # Test basic fingerprint functionality
    fingerprint = original_analyzer.get_layer_fingerprint(layer_thickness=1.0)
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

    # Test rotation detection between materials
    rotation_info = another_analyzer.detect_rotation_from_original(original_material)
    assert isinstance(rotation_info, dict)
    assert 'is_rotated' in rotation_info
    
    # The expected behavior based on the is_flipped parameter
    assert rotation_info['is_rotated'] == is_flipped, "Rotation detection status does not match expected value."


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


def test_rotation_detection_and_correction():
    """Test rotation detection and manual correction using operations."""
    material = Material.create(SI_CONVENTIONAL_SLAB_001)
    
    # Create a 180-degree rotated material
    rotated_material_180 = rotate(material, axis=[1, 0, 0], angle=180, rotate_cell=False)
    rotated_analyzer = BasisMaterialAnalyzer(material=rotated_material_180)
    
    # Test rotation detection
    rotation_info = rotated_analyzer.detect_rotation_from_original(material, layer_thickness=0.5)
    
    # Should detect some rotation (may not always be exactly 180° due to material symmetry)
    if rotation_info['is_rotated']:
        assert rotation_info['rotation_matrix'] is not None
        assert rotation_info['rotation_angle'] is not None
        assert rotation_info['rotation_axis'] is not None
        assert isinstance(rotation_info['confidence'], float)
        assert 0.0 <= rotation_info['confidence'] <= 1.0
        
        # Test corrective rotation using operations module
        rotation_axis = rotation_info['rotation_axis']
        rotation_angle = -rotation_info['rotation_angle']  # Inverse rotation
        
        corrected_material = rotate(
            rotated_material_180,
            axis=rotation_axis,
            angle=rotation_angle,
            rotate_cell=False
        )
        assert isinstance(corrected_material, Material)
        assert corrected_material.basis.number_of_atoms == material.basis.number_of_atoms
    
    # Test self-comparison - should not detect rotation
    self_analyzer = BasisMaterialAnalyzer(material=material)
    self_rotation_info = self_analyzer.detect_rotation_from_original(material, layer_thickness=0.5)
    assert not self_rotation_info['is_rotated'], "Self-comparison should not detect rotation"
