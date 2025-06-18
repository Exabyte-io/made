import difflib
import json
from typing import Any, Dict

from mat3ra.utils import assertion as assertion_utils
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


def normalize_pydantic_data(entity: Any) -> Any:
    """
    Centralized function to normalize any Pydantic models to plain Python data structures.
    Handles RootModels, regular BaseModels, lists, dicts, and primitives uniformly.
    """
    try:
        # Try to serialize using Pydantic's model_dump if available (v2)
        if hasattr(entity, "model_dump"):
            return entity.model_dump()
        # Fall back to dict() for Pydantic v1
        elif hasattr(entity, "dict"):
            return entity.dict()
        # Handle lists recursively
        elif isinstance(entity, list):
            return [normalize_pydantic_data(item) for item in entity]
        # Handle dicts recursively
        elif isinstance(entity, dict):
            return {k: normalize_pydantic_data(v) for k, v in entity.items()}
        # Return primitives as-is
        else:
            return entity
    except Exception:
        # If all else fails, return as-is
        return entity


def assert_two_entities_deep_almost_equal(entity1, entity2, rtol=1e-5, atol=1e-9):
    # Normalize both entities to plain Python data structures
    normalized_entity1 = normalize_pydantic_data(entity1)
    normalized_entity2 = normalize_pydantic_data(entity2)

    cleaned_dict_1 = prune_extra_keys(normalized_entity1, normalized_entity2)

    sorted_dict_1 = sort_dict_recursively(cleaned_dict_1)
    sorted_dict_2 = sort_dict_recursively(normalized_entity2)

    actual_data = json.loads(json.dumps(sorted_dict_1))
    expected_data = json.loads(json.dumps(sorted_dict_2))

    try:
        assertion_utils.assert_deep_almost_equal(expected_data, actual_data, rtol=rtol, atol=atol)
    except AssertionError as e:
        show_difference(expected_data, actual_data)
        raise e
