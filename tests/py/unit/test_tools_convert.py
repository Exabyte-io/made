import numpy as np
from ase import Atoms
from ase.build import bulk
from mat3ra.made.basis import Basis
from mat3ra.made.material import Material
from mat3ra.made.tools.convert import from_ase, from_poscar, from_pymatgen, to_ase, to_poscar, to_pymatgen
from mat3ra.utils import assertion as assertion_utils
from pymatgen.core.structure import Element, Lattice, Structure

from .fixtures import INTERFACE_PROPERTIES_JSON, INTERFACE_STRUCTURE

PYMATGEN_LATTICE = Lattice.from_parameters(a=3.84, b=3.84, c=3.84, alpha=120, beta=90, gamma=60)
PYMATGEN_STRUCTURE = Structure(PYMATGEN_LATTICE, ["Si", "Si"], [[0, 0, 0], [0.75, 0.5, 0.75]])

POSCAR_CONTENT = """Si2
1.0
   3.3489202364344242    0.0000000000000000    1.9335000000000004
   1.1163067454781415    3.1573922784475164    1.9335000000000004
   0.0000000000000000    0.0000000000000000    3.8670000000000000
Si
2
direct
   0.0000000000000000    0.0000000000000000    0.0000000000000000 Si
   0.2500000000000000    0.2500000000000000    0.2500000000000000 Si
"""


def test_to_pymatgen():
    material = Material.create(Material.default_config)
    structure = to_pymatgen(material)
    assert isinstance(structure, Structure)
    assert structure.lattice == Lattice.from_parameters(3.867, 3.867, 3.867, 60, 60, 60)
    assert structure.species == [Element("Si"), Element("Si")]
    assert (structure.frac_coords == [[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]]).all()


def test_from_pymatgen():
    material_data = from_pymatgen(PYMATGEN_STRUCTURE)
    assert material_data["lattice"]["a"] == 3.84
    assert material_data["lattice"]["alpha"] == 120
    interface_data = from_pymatgen(INTERFACE_STRUCTURE)
    actual_properties = interface_data["metadata"]["interface_properties"]
    assertion_utils.assert_deep_almost_equal(INTERFACE_PROPERTIES_JSON, actual_properties)


def test_to_poscar():
    material = Material.create(Material.default_config)
    poscar = to_poscar(material)
    assert poscar == POSCAR_CONTENT


def test_from_poscar():
    material_data = from_poscar(POSCAR_CONTENT)
    assert material_data["lattice"]["a"] == 3.867
    assert material_data["lattice"]["alpha"] == 60


def test_to_ase():
    material = Material.create(Material.default_config)
    labels_array = [{"id": 0, "value": 0}, {"id": 1, "value": 1}]
    material.basis = Basis.from_dict(**Material.default_config["basis"], labels=labels_array)
    ase_atoms = to_ase(material)
    assert isinstance(ase_atoms, Atoms)
    assert np.allclose(
        ase_atoms.get_cell().tolist(),
        [
            [3.3489202364344242, 0.0, 1.9335000000000004],
            [1.1163067454781415, 3.1573922784475164, 1.9335000000000004],
            [0.0, 0.0, 3.8670000000000000],
        ],
    )
    assert ase_atoms.get_chemical_symbols() == ["Si", "Si"]
    assert np.allclose(ase_atoms.get_scaled_positions().tolist(), [[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]])
    assert ase_atoms.get_tags().tolist() == [0, 1]


def test_from_ase():
    ase_atoms = bulk("Si")
    ase_atoms.set_tags([0, 1])
    material_data = from_ase(ase_atoms)
    assert material_data["lattice"]["a"] == 3.839589822
    assert material_data["lattice"]["alpha"] == 60
    assert material_data["basis"]["labels"] == [{"id": 0, "value": 0}, {"id": 1, "value": 1}]
