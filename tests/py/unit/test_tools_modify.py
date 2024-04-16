from ase.build import bulk
from mat3ra.made.tools.modify import filter_by_label


def test_filter_by_label():
    substrate = bulk("Si", cubic=True)
    film = bulk("Cu", cubic=True)
    interface_atoms = substrate + film
    interface_atoms.set_tags([1] * len(substrate) + [2] * len(film))
    film_extracted = filter_by_label(interface_atoms, 2)
    assert (film.symbols == film_extracted.symbols).all()
