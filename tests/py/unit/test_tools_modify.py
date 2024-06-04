from ase.build import bulk
from mat3ra.made.material import Material
from mat3ra.made.tools.convert import from_ase
from mat3ra.made.tools.modify import filter_by_label


def test_filter_by_label():
    substrate = bulk("Si", cubic=True)
    film = bulk("Cu", cubic=True)
    interface_atoms = substrate + film
    interface_atoms.set_tags([1] * len(substrate) + [2] * len(film))
    material_interface = Material(from_ase(interface_atoms))
    film_extracted = filter_by_label(material_interface, 2)
    film_material = Material(from_ase(film))

    # Ids of filtered elements will be missing, comparing the resulting values
    assert [el["value"] for el in film_material.basis["elements"]] == [
        el["value"] for el in film_extracted.basis["elements"]
    ]
