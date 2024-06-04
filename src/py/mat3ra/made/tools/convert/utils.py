from typing import Any

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


def extract_interface_labels_from_pymatgen(structure: type(PymatgenInterface)) -> list[dict[str, int | Any]] | None:  # type: ignore
    interface_labels = None
    if any("interface_label" in site.properties for site in structure.sites):
        interface_labels = [
            {"id": idx, "value": INTERFACE_LABELS_MAP[site.properties["interface_label"]]}
            for idx, site in enumerate(structure.sites)
        ]
    return interface_labels
