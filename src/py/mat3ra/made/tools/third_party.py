from ase import Atoms as ASEAtoms
from ase.build import add_vacuum as ase_add_vacuum
from ase.build.supercells import make_supercell as ase_make_supercell
from ase.calculators.calculator import Calculator as ASECalculator
from ase.calculators.emt import EMT as ASECalculatorEMT
from pymatgen.analysis.defects.core import Interstitial as PymatgenInterstitial
from pymatgen.analysis.defects.core import Substitution as PymatgenSubstitution
from pymatgen.analysis.defects.core import Vacancy as PymatgenVacancy
from pymatgen.analysis.local_env import VoronoiNN as PymatgenVoronoiNN
from pymatgen.core import IStructure as PymatgenIStructure
from pymatgen.core import PeriodicSite as PymatgenPeriodicSite
from pymatgen.core.interface import Interface as PymatgenInterface
from pymatgen.core.interface import label_termination as label_pymatgen_slab_termination
from pymatgen.core.structure import Lattice as PymatgenLattice
from pymatgen.core.structure import Structure as PymatgenStructure
from pymatgen.core.surface import Slab as PymatgenSlab
from pymatgen.core.surface import SlabGenerator as PymatgenSlabGenerator
from pymatgen.io.ase import AseAtomsAdaptor as PymatgenAseAtomsAdaptor
from pymatgen.io.vasp.inputs import Poscar as PymatgenPoscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer as PymatgenSpacegroupAnalyzer

# Re-exported imports to allow for both use in type hints and instantiation

__all__ = [
    "ASEAtoms",
    "ASECalculator",
    "ASECalculatorEMT",
    "PymatgenLattice",
    "PymatgenStructure",
    "PymatgenIStructure",
    "PymatgenSlab",
    "PymatgenSlabGenerator",
    "PymatgenInterface",
    "PymatgenPeriodicSite",
    "PymatgenSpacegroupAnalyzer",
    "PymatgenVacancy",
    "PymatgenSubstitution",
    "PymatgenInterstitial",
    "label_pymatgen_slab_termination",
    "ase_make_supercell",
    "ase_add_vacuum",
    "PymatgenAseAtomsAdaptor",
    "PymatgenPoscar",
    "PymatgenVoronoiNN",
]
