import json
from typing import Any, Dict, List, Union

import numpy as np
from mat3ra.utils.object import NumpyNDArrayRoundEncoder
from scipy.spatial.distance import pdist

from mat3ra.made.utils import map_array_to_array_with_id_value, get_center_of_coordinates
from .interface_parts_enum import INTERFACE_LABELS_MAP
from ..third_party import ASEAtoms, PymatgenInterface, PymatgenStructure


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


def calculate_molecule_padding_cell(coordinates: List[List[float]], padding_factor: float = 2.0) -> List[List[float]]:
    """
    Calculate values for a padded cell for a molecule based on its coordinates.
    Args:
        coordinates (Array[Array[float]]): A list of atomic coordinates.
        padding_factor (float): The factor by which to multiply the maximum distance for padding.
    Returns:
        Array[float]: A list containing the final cell latice vectors with padding applied.
    """
    positions = np.array(coordinates)
    center = get_center_of_coordinates(positions)
    shifted_positions = positions - center
    max_distance = np.max(pdist(shifted_positions)) if len(positions) >= 2 else 10.0
    padding_value = padding_factor * max_distance

    return [
        [padding_value, 0.0, 0.0],
        [0.0, padding_value, 0.0],
        [0.0, 0.0, padding_value],
    ]
