import difflib
import json
import sys
from enum import Enum
from typing import Any, Dict

from mat3ra.made.tools.third_party import PymatgenAseAtomsAdaptor
from mat3ra.made.tools.utils import unwrap
from mat3ra.utils import assertion as assertion_utils
from pymatgen.core.structure import Structure


class OSPlatform(Enum):
    """Platform enum for architecture-specific test handling."""

    DARWIN = "darwin"
    OTHER = "other"


def get_current_platform() -> OSPlatform:
    return OSPlatform.DARWIN if sys.platform == "darwin" else OSPlatform.OTHER


def get_platform_specific_value(platform_values: Any) -> Any:
    """
    Get platform-specific value for the test cases that are architecture-specific.

    Args:
        platform_values: Either a dictionary mapping TestPlatform to values, or any other value

    Returns:
        Platform-specific value if input is a platform dictionary, otherwise returns input unchanged
    """
    if isinstance(platform_values, dict) and OSPlatform.DARWIN in platform_values:
        current_platform = get_current_platform()
        return platform_values[current_platform]
    return platform_values


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

    adaptor = PymatgenAseAtomsAdaptor()
    interface_structure = adaptor.get_structure(atoms)
    interface_structure.add_site_property(
        "interface_label",
        [ATOMS_TAGS_TO_INTERFACE_STRUCTURE_LABELS[tag] for tag in interface_structure.site_properties["tags"]],
    )
    return interface_structure


def prune_extra_keys(data: Any, reference: Any) -> Any:
    if isinstance(data, dict) and isinstance(reference, dict):
        return {key: prune_extra_keys(data[key], reference[key]) for key in data if key in reference}
    elif isinstance(data, list) and isinstance(reference, list):
        return [prune_extra_keys(d_item, r_item) for d_item, r_item in zip(data, reference)]
    else:
        return data


def sort_dict_recursively(d):
    """Sort dictionary keys recursively."""
    if isinstance(d, dict):
        return {k: sort_dict_recursively(v) for k, v in sorted(d.items())}
    elif isinstance(d, list):
        return [sort_dict_recursively(x) for x in d]
    else:
        return d


def show_difference(expected_data, actual_data):
    diff = difflib.ndiff(
        json.dumps(expected_data, indent=4).splitlines(),
        json.dumps(actual_data, indent=4).splitlines(),
    )
    diff_str = "\n".join(diff)
    # Filter out the lines that are the same
    diff_str = "\n".join(line for line in diff_str.splitlines() if not line.startswith("  "))
    if diff_str:
        raise AssertionError(f"Entities differ:\n{diff_str}")


def assert_two_entities_deep_almost_equal(entity1, entity2, rtol=1e-5, atol=1e-9):
    # First unwrap any nested schema objects
    entity1 = unwrap(entity1)
    entity2 = unwrap(entity2)

    dict_1 = entity1 if isinstance(entity1, (dict, list)) else json.loads(entity1.to_json())
    dict_2 = entity2 if isinstance(entity2, (dict, list)) else json.loads(entity2.to_json())

    cleaned_dict_1 = prune_extra_keys(dict_1, dict_2)

    sorted_dict_1 = sort_dict_recursively(cleaned_dict_1)
    sorted_dict_2 = sort_dict_recursively(dict_2)

    actual_data = json.loads(json.dumps(sorted_dict_1))
    expected_data = json.loads(json.dumps(sorted_dict_2))

    try:
        assertion_utils.assert_deep_almost_equal(expected_data, actual_data, rtol=rtol, atol=atol)
    except AssertionError as e:
        show_difference(expected_data, actual_data)
        raise e
