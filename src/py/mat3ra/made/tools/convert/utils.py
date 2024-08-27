import json
from enum import Enum
from typing import Any, Dict, List, Union

from mat3ra.made.utils import map_array_to_array_with_id_value
from mat3ra.utils.object import NumpyNDArrayRoundEncoder

from ..third_party import ASEAtoms, PymatgenInterface, PymatgenStructure


class InterfacePartsEnum(str, Enum):
    SUBSTRATE = 0
    FILM = 1


INTERFACE_LABELS_MAP = {"substrate": 0, "film": 1}


def extract_labels_from_pymatgen_structure(structure: PymatgenStructure) -> List[int]:
    labels = []
    if isinstance(structure, PymatgenInterface):
        labels = list(map(lambda s: INTERFACE_LABELS_MAP[s.properties["interface_label"]], structure.sites))
    return labels


def extract_metadata_from_pymatgen_structure(structure: PymatgenStructure) -> Dict[str, Any]:
    metadata = {}
    # TODO: consider using Interface JSONSchema from ESSE when such created and adapt interface_properties accordingly.
    # Add interface properties to metadata according to pymatgen Interface as a JSON object
    if hasattr(structure, "interface_properties"):
        interface_props = structure.interface_properties
        # TODO: figure out how to round the values and stringify terminations tuple
        #  in the interface properties with Encoder
        for key, value in interface_props.items():
            if isinstance(value, tuple):
                interface_props[key] = str(value)
        metadata["interface_properties"] = json.loads(json.dumps(interface_props, cls=NumpyNDArrayRoundEncoder))

    return metadata


def extract_tags_from_ase_atoms(atoms: ASEAtoms) -> List[Union[str, int]]:
    result = []
    if "tags" in atoms.arrays:
        int_tags = [int(tag) for tag in atoms.arrays["tags"] if tag is not None]
        result = map_array_to_array_with_id_value(int_tags, remove_none=True)
    return result
