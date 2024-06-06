from typing import Any, Dict, List

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


# TODO: move to a more general location
def map_array_to_array_with_id_value(array: List[Any], remove_none: bool = False) -> List[Any]:
    full_array = [{"id": i, "value": item} for i, item in enumerate(array)]
    if remove_none:
        return list(filter(lambda x: x["value"] is not None, full_array))
    return full_array


def map_array_with_id_value_to_array(array: List[Dict[str, Any]]) -> List[Any]:
    return [item["value"] for item in array]
