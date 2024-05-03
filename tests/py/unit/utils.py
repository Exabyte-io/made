from typing import Dict

from pymatgen.core.structure import Structure
from pymatgen.io.ase import AseAtomsAdaptor

ATOMS_TAGS_TO_INTERFACE_STRUCTURE_LABELS: Dict = {1: "substrate", 2: "film"}
INTERFACE_STRUCTURE_LABELS_TO_ATOMS_TAGS: Dict = {v: k for k, v in ATOMS_TAGS_TO_INTERFACE_STRUCTURE_LABELS.items()}


def atoms_to_interface_structure(atoms) -> Structure:
    """
    Converts ASE Atoms object to pymatgen Interface object.
    Args:
        atoms (Atoms): The ASE Atoms object.
    Returns:
        Interface: The pymatgen Interface object.
    """

    adaptor = AseAtomsAdaptor()
    interface_structure = adaptor.get_structure(atoms)
    interface_structure.add_site_property(
        "interface_label",
        [ATOMS_TAGS_TO_INTERFACE_STRUCTURE_LABELS[tag] for tag in interface_structure.site_properties["tags"]],
    )
    return interface_structure
