from mat3ra.made.basis import Basis
from mat3ra.made.lattice import Lattice
from mat3ra.made.material import Material
from mat3ra.utils import assertion as assertion_utils

REFERENCE_OBJECT_1 = {"key1": "value1", "key2": "value2"}


def test_create():
    material = Material.create(Material.default_config)
    assert isinstance(material.basis, Basis)
    assert isinstance(material.lattice, Lattice)
    assert material.name == Material.default_config["name"]


def test_material_to_json():
    material = Material.create(Material.default_config)
    assertion_utils.assert_deep_almost_equal(Material.default_config, material.to_json())


def test_basis_to_json():
    material = Material.create(Material.default_config)
    basis = material.basis
    expected_basis_config = {**Material.default_config["basis"], "labels": []}
    assertion_utils.assert_deep_almost_equal(expected_basis_config, basis.to_json())


# TODO: Add test to check if basis.cell is changed when lattice of material is changed, and vice versa


def test_create_empty():
    """Test creating empty materials with different lattice parameters"""
    # Test default parameters
    material = Material.create_empty()
    assert material.basis.elements.values == []
    assert material.basis.coordinates.values == []
    assert material.lattice.type == "CUB"
    assert material.lattice.a == 1.0
    assert material.lattice.b == 1.0
    assert material.lattice.c == 1.0
    assert material.lattice.alpha == 90.0
    assert material.lattice.beta == 90.0
    assert material.lattice.gamma == 90.0
    assert material.name == "New Material"

    # Test custom parameters
    material = Material.create_empty(
        a=2.0,
        b=3.0,
        c=4.0,
        alpha=80.0,
        beta=85.0,
        gamma=95.0,
        lattice_type="TRI",
        name="Custom Empty"
    )
    assert material.basis.elements.values == []
    assert material.basis.coordinates.values == []
    assert material.lattice.type == "TRI"
    assert material.lattice.a == 2.0
    assert material.lattice.b == 3.0
    assert material.lattice.c == 4.0
    assert material.lattice.alpha == 80.0
    assert material.lattice.beta == 85.0
    assert material.lattice.gamma == 95.0
    assert material.name == "Custom Empty"

    # Test with only 'a' parameter
    material = Material.create_empty(a=2.5)
    assert material.lattice.a == 2.5
    assert material.lattice.b == 2.5
    assert material.lattice.c == 2.5

