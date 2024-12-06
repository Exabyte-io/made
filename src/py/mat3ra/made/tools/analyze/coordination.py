from typing import Dict, List

import numpy as np
from mat3ra.made.material import Material
from mat3ra.made.tools.analyze import get_nearest_neighbors_vectors
from pydantic import BaseModel, Field


class CoordinationAnalyzer(BaseModel):
    """
    Class to handle coordination number analysis for atoms in a material.
    Including vectors to nearest neighbors and identifying undercoordinated atoms.

    Args:
        cutoff (float): Cutoff radius for defining neighbors.
        coordination_threshold (int): Minimum coordination number for an atom to be considered fully coordinated.
    """

    cutoff: float = 3.0
    coordination_threshold: int = 3

    def get_coordination_numbers(self, material: Material) -> Dict[int, int]:
        """
        Calculate the coordination numbers for all atoms in the material.

        Args:
            material (Material): Material object.

        Returns:
            Dict[int, int]: A dictionary mapping atom indices to their coordination numbers.
        """
        nearest_neighbors = get_nearest_neighbors_vectors(material=material, cutoff=self.cutoff)
        coordination_numbers = {idx: len(vectors) for idx, vectors in enumerate(nearest_neighbors.values)}
        return coordination_numbers

    def get_undercoordinated_atom_indices(self, material: Material) -> List[int]:
        """
        Identify undercoordinated atoms based on the coordination threshold (inclusive).

        Args:
            material (Material): Material object.

        Returns:
            List[int]: List of indices of undercoordinated atoms.
        """
        coordination_numbers = self.get_coordination_numbers(material)
        return [idx for idx, number in coordination_numbers.items() if number <= self.coordination_threshold]

    def get_unique_coordination_numbers(self, material: Material) -> List[int]:
        """
        Get the unique coordination numbers for all atoms in the material.

        Args:
            material (Material): Material object.

        Returns:
            Set[int]: A set of unique coordination numbers present in the material.
        """
        coordination_numbers = self.get_coordination_numbers(material)
        return sorted(list(set(coordination_numbers.values())))

    angle_tolerance: float = Field(0.1, description="Tolerance for comparing angles between bond vectors.")
    max_bonds_to_passivate: int = Field(
        2, description="Maximum number of bonds to passivate for each undercoordinated atom."
    )

    @staticmethod
    def are_bonds_templates_similar(template1: np.ndarray, template2: np.ndarray, tolerance: float = 0.1) -> bool:
        """
        Check if two bond templates are similar.

        Args:
            template1 (np.ndarray): First template of bond vectors.
            template2 (np.ndarray): Second template of bond vectors.
            tolerance (float): Angle tolerance for comparison.

        Returns:
            bool: True if the templates are similar, False otherwise.
        """
        if len(template1) != len(template2):
            return False

        dot_matrix = np.dot(template1, template2.T)
        norms1 = np.linalg.norm(template1, axis=1)
        norms2 = np.linalg.norm(template2, axis=1)
        cosine_matrix = dot_matrix / np.outer(norms1, norms2)
        angles_matrix = np.arccos(np.clip(cosine_matrix, -1.0, 1.0))

        unmatched = list(range(len(template2)))
        for angle_row in angles_matrix:
            matches = np.where(angle_row < tolerance)[0]
            if len(matches) == 0:
                return False
            unmatched.remove(matches.tolist()[0])

        return True

    def find_template_vectors(
        self, atom_vectors: List[List[np.ndarray]], atom_elements: List[str]
    ) -> Dict[str, List[np.ndarray]]:
        """
        Find unique bond templates for each element type.

        Args:
            atom_vectors (List[List[np.ndarray]]): List of bond vectors for each atom.
            atom_elements (List[str]): List of chemical elements for each atom.

        Returns:
            Dict[str, List[np.ndarray]]: Dictionary mapping element types to unique bond templates.
        """
        element_templates = {}

        for element in set(atom_elements):
            element_indices = [i for i, e in enumerate(atom_elements) if e == element]
            element_vector_lists = [np.array(atom_vectors[i]) for i in element_indices]

            if not element_vector_lists:
                continue

            max_coord = max(len(vectors) for vectors in element_vector_lists)
            max_coord_vectors = [v for v in element_vector_lists if len(v) == max_coord]

            unique_templates: List[np.ndarray] = []
            for template in max_coord_vectors:
                if not any(
                    self.are_bonds_templates_similar(template, existing, self.angle_tolerance)
                    for existing in unique_templates
                ):
                    unique_templates.append(template)

            element_templates[element] = unique_templates

        return element_templates

    def reconstruct_missing_bonds(
        self,
        nearest_neighbor_vectors: List[List[np.ndarray]],
        chemical_elements: List[str],
        templates: Dict[str, List[np.ndarray]],
    ) -> Dict[int, List[List[float]]]:
        """
        Reconstruct missing bonds for undercoordinated atoms.

        Args:
            nearest_neighbor_vectors (List[List[np.ndarray]]): List of bond vectors for each atom.
            chemical_elements (List[str]): List of chemical elements for each atom.
            templates (Dict[str, List[np.ndarray]]): Dictionary of bond templates for each element.

        Returns:
            Dict[int, List[List[float]]]: Dictionary mapping atom indices to reconstructed bond vectors.
        """
        missing_bonds = {}

        for idx, (vectors, element) in enumerate(zip(nearest_neighbor_vectors, chemical_elements)):
            if element not in templates:
                continue

            existing_vectors = np.array(vectors) if vectors else np.empty((0, 3))
            max_coordination = len(templates[element][0])

            if len(existing_vectors) >= max_coordination:
                continue

            best_missing = None
            best_match_count = -1

            for template in templates[element]:
                if existing_vectors.size == 0:
                    match_count = 0
                else:
                    dot_matrix = np.dot(template, existing_vectors.T)
                    cosine_matrix = dot_matrix / (
                        np.linalg.norm(template, axis=1)[:, None] * np.linalg.norm(existing_vectors, axis=1)
                    )
                    angles_matrix = np.arccos(np.clip(cosine_matrix, -1.0, 1.0))

                    matches = np.any(angles_matrix < self.angle_tolerance, axis=1)
                    match_count = np.sum(matches)

                missing = template[~matches] if existing_vectors.size != 0 else template

                if match_count > best_match_count:
                    best_match_count = match_count
                    best_missing = missing

            if best_missing is not None:
                num_bonds_to_add = min(
                    len(best_missing),
                    self.max_bonds_to_passivate,
                    max_coordination - len(existing_vectors),
                )
                missing_bonds[idx] = best_missing[:num_bonds_to_add].tolist()

        return missing_bonds
