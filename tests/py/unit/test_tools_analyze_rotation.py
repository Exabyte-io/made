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
    assert hasattr(rotation_info, "angle")
    assert hasattr(rotation_info, "axis")
    assert hasattr(rotation_info, "confidence")

    assert rotation_info.is_rotated == is_flipped, "Rotation detection status does not match expected value."


@pytest.mark.parametrize(
    "rotation_angle, rotation_axis",
    [
        (0, [1, 0, 0]),
        (90, [1, 0, 0]),
        (180, [1, 0, 0]),
        (-90, [1, 0, 0]),
        (90, [0, 0, 1]),
    ],
)
def test_rotation_detection_and_correction(rotation_angle, rotation_axis):
    """Test that rotation detection runs without errors for various rotations."""
    material = Material.create(SI_CONVENTIONAL_SLAB_001)

    rotated_material = rotate(material, axis=rotation_axis, angle=rotation_angle, rotate_cell=False)

    rotation_analyzer = MaterialRotationAnalyzer(material=rotated_material)

    rotation_info = rotation_analyzer.detect_rotation_from_original(material, layer_thickness=0.5)

    assert rotation_info.rotation_matrix is not None
    assert rotation_info.angle is not None
    assert rotation_info.axis is not None
    assert isinstance(rotation_info.confidence, float)
    assert 0.0 <= rotation_info.confidence <= 1.0
    assert isinstance(rotation_info.is_rotated, bool)

    if rotation_info.is_rotated:
        detected_axis = rotation_info.axis
        detected_angle = -rotation_info.angle

        corrected_material = rotate(rotated_material, axis=detected_axis, angle=detected_angle, rotate_cell=False)
        assert isinstance(corrected_material, Material)
        assert corrected_material.basis.number_of_atoms == material.basis.number_of_atoms
