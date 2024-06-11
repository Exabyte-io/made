from ase import Atoms

from pymatgen.core.interface import Interface, label_termination
from pymatgen.core.structure import Lattice, Structure
from pymatgen.core.surface import Slab

PymatgenLattice = Lattice
PymatgenStructure = Structure
PymatgenSlab = Slab
PymatgenInterface = Interface
ASEAtoms = Atoms
label_pymatgen_slab_termination = label_termination

INTERFACE_LABELS_MAP = {"substrate": 0, "film": 1}
