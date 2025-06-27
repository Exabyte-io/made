import difflib
import json
from typing import Any, Dict, List, Tuple

import numpy as np
from mat3ra.made.material import Material
from mat3ra.made.tools.analyze.material import MaterialWithCrystalSites
from mat3ra.made.tools.analyze.other import get_atomic_coordinates_extremum
from mat3ra.made.tools.utils import unwrap
from mat3ra.utils import assertion as assertion_utils
from pydantic import BaseModel, Field
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


class MaterialStructureComparator(BaseModel):
    """
    Compare structural equivalence of two materials using neighbor vector analysis.

    This class leverages the existing Material infrastructure and analyzes the local
    atomic environment by computing vectors to nearest neighbors, then compares these
    vector patterns to determine if structures are equivalent under various transformations.

    The class uses composition with MaterialWithCrystalSites rather than inheritance,
    as it's a tool for analyzing materials rather than being a material itself.
    """

    neighbor_cutoff: float = Field(default=5.0, description="Maximum distance to consider for neighbors (Angstroms)")
    max_neighbors: int = Field(default=12, description="Maximum number of neighbors to analyze per atom")
    vector_tolerance: float = Field(default=0.1, description="Relative tolerance for comparing vector lengths")
    angle_tolerance: float = Field(default=5.0, description="Tolerance for comparing angles between vectors (degrees)")
    allow_c_axis_scaling: bool = Field(
        default=True, description="Whether to allow different c-axis lengths with vacuum"
    )
    min_vacuum_for_c_flexibility: float = Field(default=10.0, description="Minimum vacuum to allow c-axis flexibility")

    class Config:
        arbitrary_types_allowed = True

    def ensure_material(self, entity) -> Material:
        if isinstance(entity, Material):
            return entity

        entity = unwrap(entity)
        if isinstance(entity, dict):
            return Material.create(entity)
        elif hasattr(entity, "to_json"):
            return Material.create(json.loads(entity.to_json()))
        else:
            raise ValueError(f"Cannot convert entity to Material: {type(entity)}")

    def create_material_with_sites(self, material: Material) -> MaterialWithCrystalSites:
        return MaterialWithCrystalSites.from_material(material)

    def get_neighbor_patterns(self, material: Material) -> List[Dict]:
        material_with_sites = self.create_material_with_sites(material)
        neighbor_vectors = material_with_sites.get_neighbors_vectors_for_all_sites(
            cutoff=self.neighbor_cutoff, max_number_of_neighbors=self.max_neighbors
        )

        patterns = []
        for i, site_vectors in enumerate(neighbor_vectors.values):
            if len(site_vectors) < 2:
                continue

            element = material.basis.elements.get_element_value_by_index(i)
            distances = [np.linalg.norm(vec) for vec in site_vectors]

            angles = []
            for j in range(len(site_vectors)):
                for k in range(j + 1, len(site_vectors)):
                    v1, v2 = site_vectors[j], site_vectors[k]
                    cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
                    cos_angle = np.clip(cos_angle, -1, 1)
                    angle = np.degrees(np.arccos(cos_angle))
                    angles.append(angle)

            neighbor_elements = []
            site_coord = material.basis.coordinates.get_element_value_by_index(i)

            for vec in site_vectors:
                neighbor_coord = np.array(site_coord) + vec

                min_dist = float("inf")
                closest_element = None

                for j, coord in enumerate(material.basis.coordinates.values):
                    dist = np.linalg.norm(np.array(coord) - neighbor_coord)
                    if dist < min_dist:
                        min_dist = dist
                        closest_element = material.basis.elements.get_element_value_by_index(j)

                if closest_element:
                    neighbor_elements.append(closest_element)

            pattern = {
                "element": element,
                "distances": sorted(distances),
                "angles": sorted(angles),
                "neighbor_elements": sorted(neighbor_elements),
            }

            patterns.append(pattern)

        return patterns

    def compare_patterns(self, patterns1: List[Dict], patterns2: List[Dict]) -> bool:
        if len(patterns1) != len(patterns2):
            return False

        def pattern_key(p):
            return (p["element"], tuple(p["distances"][:3]))

        patterns1_sorted = sorted(patterns1, key=pattern_key)
        patterns2_sorted = sorted(patterns2, key=pattern_key)

        for p1, p2 in zip(patterns1_sorted, patterns2_sorted):
            if not self.patterns_match(p1, p2):
                return False

        return True

    def patterns_match(self, p1: Dict, p2: Dict) -> bool:
        if p1["element"] != p2["element"]:
            return False

        if p1["neighbor_elements"] != p2["neighbor_elements"]:
            return False

        distances1, distances2 = p1["distances"], p2["distances"]
        if len(distances1) != len(distances2):
            return False

        for i, (d1, d2) in enumerate(zip(distances1, distances2)):
            rel_diff = abs(d1 - d2) / max(d1, d2)
            if rel_diff > self.vector_tolerance:
                return False

        angles1, angles2 = p1["angles"], p2["angles"]
        if len(angles1) != len(angles2):
            return False

        for i, (a1, a2) in enumerate(zip(angles1, angles2)):
            abs_diff = abs(a1 - a2)
            if abs_diff > self.angle_tolerance:
                return False

        return True

    def has_sufficient_vacuum(self, material: Material) -> bool:
        min_z = get_atomic_coordinates_extremum(material, extremum="min", axis="z", use_cartesian_coordinates=True)
        max_z = get_atomic_coordinates_extremum(material, extremum="max", axis="z", use_cartesian_coordinates=True)

        z_range = max_z - min_z

        c_length = material.lattice.c
        vacuum = c_length - z_range

        return vacuum >= self.min_vacuum_for_c_flexibility

    def _basic_compatibility_checks(self, material1: Material, material2: Material) -> bool:
        if material1.basis.number_of_atoms != material2.basis.number_of_atoms:
            return False

        elements1 = sorted(material1.basis.elements.values)
        elements2 = sorted(material2.basis.elements.values)
        if elements1 != elements2:
            return False

        return True

    def _handle_c_axis_scaling(self, material1: Material, material2: Material) -> Material:
        if not self.allow_c_axis_scaling:
            return material1

        has_vacuum1 = self.has_sufficient_vacuum(material1)
        has_vacuum2 = self.has_sufficient_vacuum(material2)

        if has_vacuum1 or has_vacuum2:
            c_ratio = material2.lattice.c / material1.lattice.c
            material1_scaled = material1.clone()
            material1_scaled.to_cartesian()

            new_coordinates = []
            for coord in material1_scaled.basis.coordinates.values:
                new_coord = [coord[0], coord[1], coord[2] * c_ratio]
                new_coordinates.append(new_coord)
            material1_scaled.basis.coordinates.values = new_coordinates
            return material1_scaled

        return material1

    def _handle_gamma_equivalent_structures(self, material1: Material, material2: Material):
        gamma1 = material1.lattice.gamma
        gamma2 = material2.lattice.gamma
        gamma_diff = abs(gamma1 - gamma2)
        is_gamma_equivalent = (
            abs(gamma_diff - 60) <= self.angle_tolerance or abs(gamma_diff - 120) <= self.angle_tolerance
        )
        if is_gamma_equivalent:
            self.vector_tolerance = 0.5
            self.angle_tolerance = 25.0

    def _get_cartesian_materials(self, material1: Material, material2: Material) -> Tuple[Material, Material]:
        material1_prepared = material1.clone()
        material2_prepared = material2.clone()
        material1_prepared.to_cartesian()
        material2_prepared.to_cartesian()
        return material1_prepared, material2_prepared

    def _compare_neighbor_patterns(self, material1: Material, material2: Material) -> bool:
        patterns1 = self.get_neighbor_patterns(material1)
        patterns2 = self.get_neighbor_patterns(material2)

        return self.compare_patterns(patterns1, patterns2)

    def compare(self, entity1, entity2) -> bool:
        """
        Compare two materials for structural equivalence using existing Material infrastructure.

        Args:
            entity1, entity2: Material objects or dictionaries to compare

        Returns:
            bool: True if structures are equivalent, False otherwise
        """
        try:
            material1 = self.ensure_material(entity1)
            material2 = self.ensure_material(entity2)

            if not self._basic_compatibility_checks(material1, material2):
                return False

            material1 = self._handle_c_axis_scaling(material1, material2)
            self._handle_gamma_equivalent_structures(material1, material2)
            material1, material2 = self._get_cartesian_materials(material1, material2)

            return self._compare_neighbor_patterns(material1, material2)

        except Exception:
            return False


def assert_slab_structures_almost_equal(
    entity1,
    entity2,
    angle_tolerance=5.0,
    c_axis_flexible=True,
    min_vacuum_for_c_flexibility=10.0,
    neighbor_cutoff=5.0,
    max_neighbors=12,
    vector_tolerance=0.05,
):
    """
    Compare slab-like structures using MaterialStructureComparator.

    Args:
        entity1, entity2: Materials or material configs to compare
        angle_tolerance: Absolute tolerance for lattice angles (degrees)
        c_axis_flexible: Whether c-axis can differ if sufficient vacuum
        min_vacuum_for_c_flexibility: Minimum vacuum size to allow c-axis flexibility
        neighbor_cutoff: Maximum distance for neighbor search
        max_neighbors: Maximum number of neighbors per atom
        vector_tolerance: Relative tolerance for bond distances (5%)

    Raises:
        AssertionError: If structures are not equivalent
    """
    comparator = MaterialStructureComparator(
        neighbor_cutoff=neighbor_cutoff,
        max_neighbors=max_neighbors,
        vector_tolerance=vector_tolerance,
        angle_tolerance=angle_tolerance,
        allow_c_axis_scaling=c_axis_flexible,
        min_vacuum_for_c_flexibility=min_vacuum_for_c_flexibility,
    )

    result = comparator.compare(entity1, entity2)
    if not result:
        raise AssertionError("Structures are not equivalent")
