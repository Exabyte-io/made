from mat3ra.made.cell import Cell

VECTORS = [
    [1.0, 0.0, 0.0],
    [0.0, 2.0, 0.0],
    [0.0, 0.0, 3.0],
]


def test_cell_creation():
    cell = Cell.from_vectors_array(vectors=VECTORS)
    assert cell.a.value == VECTORS[0]
    assert cell.b.value == VECTORS[1]
    assert cell.c.value == VECTORS[2]
    assert cell.volume == 6.0
    assert cell.vector_arrays == VECTORS


def test_cell_creation_default():
    cell = Cell()
    assert cell.a.value == Cell.__default_vectors__[0]
    assert cell.b.value == Cell.__default_vectors__[1]
    assert cell.c.value == Cell.__default_vectors__[2]
    assert cell.volume == 1.0
    assert cell.vector_arrays == Cell.__default_vectors__
