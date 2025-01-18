from mat3ra.made.lattice import Lattice


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
