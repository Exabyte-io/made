from typing import Dict

from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.core.interface import Interface

tags_to_labels: Dict = {1: "substrate", 2: "film"}
labels_to_tags: Dict = {v: k for k, v in tags_to_labels.items()}


def atoms_to_interface_structure(atoms):
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
        "interface_label", [tags_to_labels[tag] for tag in interface_structure.site_properties["tags"]]
    )
    return interface_structure
