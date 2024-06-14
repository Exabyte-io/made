from typing import List, Dict

import numpy as np
from mat3ra.made.material import Material


def merge_materials(materials: List[Material], distance_tolerance: float = 0.01) -> Material:
    """
    Merge multiple materials into a single material, replacing colliding atoms with the latest material's atoms.
    """
    if not materials:
        return None

    merged_items: Dict = {}

    for material_index, material in enumerate(materials):
        for elem, coord in zip(material.basis["elements"], material.basis["coordinates"]):
            np_coord = np.array(coord["value"])
            item_id = elem["id"]
            duplicate_found = False

            for existing_id, existing in list(merged_items.items()):
                existing_coord = np.array(existing["coordinate"]["value"])
                if np.linalg.norm(np_coord - existing_coord) < distance_tolerance:
                    if material_index > merged_items[existing_id].get("material_index", -1):
                        merged_items[existing_id] = {
                            "element": elem,
                            "coordinate": coord,
                            "material_index": material_index,
                        }
                        duplicate_found = True
                        break

            if not duplicate_found:
                merged_items[item_id] = {"element": elem, "coordinate": coord, "material_index": material_index}

    final_elements = [item["element"] for item in merged_items.values()]
    final_coordinates = [item["coordinate"] for item in merged_items.values()]

    final_labels = []
    if any("labels" in mat.basis for mat in materials):
        for material in materials:
            if "labels" in material.basis:
                for label in material.basis["labels"]:
                    if label["id"] in merged_items:
                        final_labels.append(label)

    merged_material_config = {
        "name": "Merged Material",
        "lattice": materials[0].lattice,
        "basis": {"elements": final_elements, "coordinates": final_coordinates, "labels": final_labels},
    }
    return Material(merged_material_config)
