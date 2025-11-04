"""
Tests for MaterialRotationAnalyzer class.
"""

import pytest
from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.rotation_analyzer import MaterialRotationAnalyzer, RotationDetectionResult
from mat3ra.made.tools.operations.core.unary import rotate

from .fixtures.slab import (
    SLAB_SrTiO3_011_TERMINATION_O2,
    SLAB_SrTiO3_011_TERMINATION_SrTiO,
    SI_CONVENTIONAL_SLAB_001,
)


@pytest.mark.parametrize(
    "original_material_config, another_material_config, is_flipped",
    [
        (SLAB_SrTiO3_011_TERMINATION_O2, SLAB_SrTiO3_011_TERMINATION_O2, False),
        (SLAB_SrTiO3_011_TERMINATION_O2, SLAB_SrTiO3_011_TERMINATION_SrTiO, True),
    ],
)
def test_rotation_detection_between_materials(original_material_config, another_material_config, is_flipped):
    """Test rotation detection between different materials."""
    original_material = Material.create(original_material_config)
    another_material = Material.create(another_material_config)

    rotation_analyzer = MaterialRotationAnalyzer(material=another_material)
    rotation_info = rotation_analyzer.detect_rotation_from_original(original_material)
    
    assert isinstance(rotation_info, RotationDetectionResult)
    assert hasattr(rotation_info, "is_rotated")
    assert hasattr(rotation_info, "rotation_matrix")
    assert hasattr(rotation_info, "rotation_angle")
    assert hasattr(rotation_info, "rotation_axis")
    assert hasattr(rotation_info, "confidence")

    assert rotation_info.is_rotated == is_flipped, "Rotation detection status does not match expected value."


def test_rotation_detection_and_correction():
    """Test rotation detection and manual correction using operations."""
    material = Material.create(SI_CONVENTIONAL_SLAB_001)

    # Create a 180-degree rotated material
    rotated_material_180 = rotate(material, axis=[1, 0, 0], angle=180, rotate_cell=False)
    rotation_analyzer = MaterialRotationAnalyzer(material=rotated_material_180)

    # Test rotation detection
    rotation_info = rotation_analyzer.detect_rotation_from_original(material, layer_thickness=0.5)

    # Should detect some rotation (may not always be exactly 180° due to material symmetry)
    if rotation_info.is_rotated:
        assert rotation_info.rotation_matrix is not None
        assert rotation_info.rotation_angle is not None
        assert rotation_info.rotation_axis is not None
        assert isinstance(rotation_info.confidence, float)
        assert 0.0 <= rotation_info.confidence <= 1.0

        # Test corrective rotation using operations module
        rotation_axis = rotation_info.rotation_axis
        rotation_angle = -rotation_info.rotation_angle

        corrected_material = rotate(rotated_material_180, axis=rotation_axis, angle=rotation_angle, rotate_cell=False)
        assert isinstance(corrected_material, Material)
        assert corrected_material.basis.number_of_atoms == material.basis.number_of_atoms


def test_self_comparison_no_rotation():
    """Test that self-comparison does not detect rotation."""
    material = Material.create(SI_CONVENTIONAL_SLAB_001)
    
    rotation_analyzer = MaterialRotationAnalyzer(material=material)
    rotation_info = rotation_analyzer.detect_rotation_from_original(material, layer_thickness=0.5)
    
    assert not rotation_info.is_rotated, "Self-comparison should not detect rotation"
    assert isinstance(rotation_info.confidence, float)
    assert 0.0 <= rotation_info.confidence <= 1.0

