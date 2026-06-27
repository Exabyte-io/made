import copy
import difflib
import json
import sys
from enum import Enum
from typing import Any, Dict

import numpy as np
from mat3ra.made.basis import Basis
from mat3ra.made.tools.third_party import PymatgenAseAtomsAdaptor
from mat3ra.made.tools.utils import unwrap
from mat3ra.utils import assertion as assertion_utils
from pymatgen.core.structure import Structure

UPDATED_COORDINATE_TOLERANCE = 1e-6


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


def to_basis_obj(b_dict: dict) -> Basis:
    """
    Ensure a basis dictionary has sequential IDs and convert it to a Basis object.
    """
    b_dict = copy.deepcopy(b_dict)
    for key in ["elements", "coordinates"]:
        if key in b_dict and isinstance(b_dict[key], list):
            for idx, item in enumerate(b_dict[key]):
                if isinstance(item, dict) and "id" not in item:
                    item["id"] = idx
    return Basis(**b_dict)


def get_basis_arrays(b: Basis) -> tuple:
    n = b.number_of_atoms
    lbls = b.labels.values if (b.labels and len(b.labels.values) == n) else [None] * n
    consts = b.constraints.values if (b.constraints and len(b.constraints.values) == n) else [None] * n
    return b.elements.values, b.coordinates.values, lbls, consts


def are_bases_almost_equal(b1: Basis, b2: Basis, atol: float) -> bool:
    """
    Compare two Basis objects for equality with permuted ordering and coordinate tolerance.
    """
    if b1.number_of_atoms != b2.number_of_atoms:
        return False
    n = b1.number_of_atoms
    self_elems, self_coords, self_lbls, self_consts = get_basis_arrays(b1)
    other_elems, other_coords, other_lbls, other_consts = get_basis_arrays(b2)

    other_indices_left = set(range(n))
    for i in range(n):
        found = False
        elem, coord, lbl, const = self_elems[i], self_coords[i], self_lbls[i], self_consts[i]
        for j in list(other_indices_left):
            if other_elems[j] == elem and other_lbls[j] == lbl and other_consts[j] == const:
                if np.allclose(other_coords[j], coord, atol=atol):
                    other_indices_left.remove(j)
                    found = True
                    break
        if not found:
            return False
    return True


def print_bases_diff(b1: Basis, b2: Basis):
    print("\n=== BASE ALMOST EQUAL FAILED ===")
    print("b1 elements:", b1.elements.values)
    print("b1 coords:", [c.tolist() if hasattr(c, "tolist") else c for c in b1.coordinates.values])
    print("b2 elements:", b2.elements.values)
    print("b2 coords:", [c.tolist() if hasattr(c, "tolist") else c for c in b2.coordinates.values])


def assert_interfaces_almost_equal(interface1: Any, interface2: Any) -> None:
    """
    Assert that two interface structures are almost equal, checking bases permutation-invariantly.
    """
    obj1 = unwrap(interface1)
    obj2 = unwrap(interface2)
    b1 = to_basis_obj(obj1.basis.to_dict())
    b2 = to_basis_obj(obj2.basis.to_dict() if hasattr(obj2, "basis") else obj2["basis"])
    if b1.units != b2.units:
        b2 = b2.clone()
        b2.to_crystal() if b1.is_in_crystal_units else b2.to_cartesian()
    atol = UPDATED_COORDINATE_TOLERANCE
    if not are_bases_almost_equal(b1, b2, atol=atol):
        print_bases_diff(b1, b2)
        assert False, "Bases are not almost equal"
    dict_1 = copy.deepcopy(json.loads(obj1.to_json()))
    dict_2 = copy.deepcopy(obj2 if isinstance(obj2, dict) else json.loads(obj2.to_json()))
    for d in [dict_1, dict_2]:
        if "basis" in d:
            del d["basis"]
        if "name" in d:
            del d["name"]
    assert_two_entities_deep_almost_equal(dict_1, dict_2, atol=atol)
