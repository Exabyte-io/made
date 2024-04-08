from mat3ra.made.material import Material
from mat3ra.made.tools.material import to_pymatgen
from pymatgen.core.structure import Element, Lattice, Structure


def test_to_pymatgen():
    material = Material.create(Material.default_config)
    structure = to_pymatgen(material)
    assert isinstance(structure, Structure)
    assert structure.lattice == Lattice.from_parameters(3.867, 3.867, 3.867, 60, 60, 60)
    assert structure.species == [Element("Si"), Element("Si")]
    assert (structure.frac_coords == [[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]]).all()
