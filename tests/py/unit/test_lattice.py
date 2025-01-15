from mat3ra.made.lattice import Lattice
from mat3ra.made.material import Material

def test_create_empty_custom_lattice():
    """Test custom lattice parameters when creating empty material"""
    material = Material.create_empty(
        a=2.0, b=3.0, c=4.0, alpha=80.0, beta=85.0, gamma=95.0, lattice_type="TRI"
    )
    assert material.lattice.type == "TRI"
    assert material.lattice.a == 2.0
    assert material.lattice.b == 3.0
    assert material.lattice.c == 4.0
    assert material.lattice.alpha == 80.0
    assert material.lattice.beta == 85.0
    assert material.lattice.gamma == 95.0


def test_lattice_vectors_access():
    lattice = Lattice(a=2.0, b=3.0, c=4.0)

    # Test individual vector access
    assert isinstance(lattice.vectors.a, list)
    assert isinstance(lattice.vectors.b, list)
    assert isinstance(lattice.vectors.c, list)

    # Test vector arrays access
    arrays = lattice.vector_arrays
    assert len(arrays) == 3
    assert arrays[0] == lattice.vectors.a
    assert arrays[1] == lattice.vectors.b
    assert arrays[2] == lattice.vectors.c
