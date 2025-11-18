from mat3ra.esse.models.core.reusable.axis_enum import AxisEnum
from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.basis import BasisMaterialAnalyzer
from mat3ra.made.tools.analyze.fingerprint import LayeredFingerprintAlongAxis
from mat3ra.made.tools.operations.core.unary import supercell

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

    non_empty_layers = fingerprint.non_empty_layers
    assert isinstance(non_empty_layers, list)
    assert all(isinstance(layer, type(fingerprint.layers[0])) for layer in non_empty_layers)

    element_sequence = fingerprint.element_sequence
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


def test_fingerprint_similarity_scores():
    """Test get_similarity_score and get_similarity_score_ignore_periodicity methods with real materials."""
    slab = Material.create(SLAB_SrTiO3_011_TERMINATION_O2)
    supercell_matrix = [[2, 0, 0], [0, 2, 0], [0, 0, 1]]
    supercell_material = supercell(slab, supercell_matrix)

    analyzer_slab = BasisMaterialAnalyzer(material=slab)
    analyzer_supercell = BasisMaterialAnalyzer(material=supercell_material)

    layer_thickness = 1.0

    # Test along y-axis where 2x2x1 supercell creates 2x repetition
    fp_slab_y = analyzer_slab.get_layer_fingerprint(layer_thickness=layer_thickness, axis=AxisEnum.y)
    fp_supercell_y = analyzer_supercell.get_layer_fingerprint(layer_thickness=layer_thickness, axis=AxisEnum.y)

    # Test get_similarity_score with different number of layers (should return 0.0)
    assert fp_slab_y.get_similarity_score(fp_supercell_y) == 0.0

    # Test get_similarity_score_ignore_periodicity - supercell should be a periodic repetition
    # The supercell has 2x repetition in y, so along y-axis it should detect periodicity
    # Note: Score may not be exactly 1.0 due to layer boundary alignment, but should be > 0
    score_periodic = fp_slab_y.get_similarity_score_ignore_periodicity(fp_supercell_y)
    assert score_periodic > 0.0
    assert len(fp_supercell_y.layers) % len(fp_slab_y.layers) == 0

    # Test bidirectional
    assert fp_supercell_y.get_similarity_score_ignore_periodicity(fp_slab_y) > 0.0

    # Test along z-axis where layers should match (no repetition in z for 2x2x1)
    fp_slab_z = analyzer_slab.get_layer_fingerprint(layer_thickness=layer_thickness, axis=AxisEnum.z)
    fp_supercell_z = analyzer_supercell.get_layer_fingerprint(layer_thickness=layer_thickness, axis=AxisEnum.z)

    # Along z-axis, both should have same number of layers and match
    assert len(fp_slab_z.layers) == len(fp_supercell_z.layers)
    assert fp_slab_z.get_similarity_score(fp_supercell_z) == 1.0
