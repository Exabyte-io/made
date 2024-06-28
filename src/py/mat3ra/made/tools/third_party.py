from ase import Atoms as ASEAtoms
from pymatgen.core import IStructure as PymatgenIStructure
from pymatgen.core import PeriodicSite as PymatgenPeriodicSite
from pymatgen.core.interface import Interface as PymatgenInterface
from pymatgen.core.interface import label_termination
from pymatgen.core.structure import Lattice as PymatgenLattice
from pymatgen.core.structure import Structure as PymatgenStructure
from pymatgen.core.surface import Slab as PymatgenSlab
from pymatgen.core.surface import SlabGenerator as PymatgenSlabGenerator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer as PymatgenSpacegroupAnalyzer

# Re-exported imports to allow for both use in type hints and instantiation
label_pymatgen_slab_termination = label_termination


__all__ = [
    "ASEAtoms",
    "PymatgenLattice",
    "PymatgenStructure",
    "PymatgenIStructure",
    "PymatgenSlab",
    "PymatgenSlabGenerator",
    "PymatgenInterface",
    "PymatgenPeriodicSite",
    "PymatgenSpacegroupAnalyzer",
    "label_pymatgen_slab_termination",
]
