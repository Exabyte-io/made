from enum import Enum
from typing import Dict, List

import numpy as np


class BondDirectionsTemplatesEnum(list, Enum):
    OCTAHEDRAL = [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]]
    TETRAHEDRAL = [[1, 1, 1], [-1, -1, 1], [1, -1, -1], [-1, 1, -1]]
    PLANAR = [[1, 0, 0], [-1, 2, 0], [-1, -2, 0]]
    LINEAR = [[1, 0, 0], [-1, 0, 0]]


class BondDirections(np.ndarray):
    def __eq__(self, other):
        if not isinstance(other, BondDirections):
            return False
        return self.are_bond_directions_similar(self, other)

    @staticmethod
    def are_bond_directions_similar(directions1: np.ndarray, directions2: np.ndarray, tolerance: float = 0.1) -> bool:
        """
        Check if two bond templates are similar.

        Args:
            directions1 (np.ndarray): First template of bond vectors.
            directions2 (np.ndarray): Second template of bond vectors.
            tolerance (float): Angle tolerance for comparison.

        Returns:
            bool: True if the templates are similar, False otherwise.
        """
        if len(directions1) != len(directions2):
            return False

        dot_matrix = np.dot(directions1, directions2.T)
        norms1 = np.linalg.norm(directions1, axis=1)
        norms2 = np.linalg.norm(directions2, axis=1)
        cosine_matrix = dot_matrix / np.outer(norms1, norms2)
        angles_matrix = np.arccos(np.clip(cosine_matrix, -1.0, 1.0))

        unmatched = list(range(len(directions2)))
        for angle_row in angles_matrix:
            matches = np.where(angle_row < tolerance)[0]
            if len(matches) == 0:
                return False
            unmatched.remove(matches.tolist()[0])

        return True

    @staticmethod
    def find_missing_directions(
        vectors: List[np.ndarray],
        # TODO: remove element
        element: str,
        # TODO: single template
        templates: List[np.ndarray],
        angle_tolerance: float = 0.1,
        max_bonds_to_passivate: int = 1,
    ) -> List[List[float]]:
        """
        Reconstruct missing bonds for a single element.

        Args:
            vectors (List[np.ndarray]): List of bond vectors for the atom.
            element (str): Chemical element of the atom.
            templates (Dict[str, List[np.ndarray]]): Dictionary of bond templates for each element.
            angle_tolerance (float): Tolerance for comparing angles between bond vectors.
            max_bonds_to_passivate (int): Maximum number of bonds to passivate for the atom.

        Returns:
            List[List[float]]: List of reconstructed bond vectors.
        """
        if element not in templates:
            return []

        existing_vectors = np.array(vectors) if vectors else np.empty((0, 3))
        max_coordination_number = len(templates[element][0])

        if len(existing_vectors) >= max_coordination_number:
            return []

        best_missing = None
        best_match_count = -1

        for template in templates[element]:
            if existing_vectors.size == 0:
                match_count = 0
            else:
                # TODO: optimize
                dot_matrix = np.dot(template, existing_vectors.T)
                cosine_matrix = dot_matrix / (
                    np.linalg.norm(template, axis=1)[:, None] * np.linalg.norm(existing_vectors, axis=1)
                )
                angles_matrix = np.arccos(np.clip(cosine_matrix, -1.0, 1.0))

                matches = np.any(angles_matrix < angle_tolerance, axis=1)
                match_count = np.sum(matches)

            missing = template[~matches] if existing_vectors.size != 0 else template

            if match_count > best_match_count:
                best_match_count = match_count
                best_missing = missing

        if best_missing is not None:
            num_bonds_to_add = min(
                len(best_missing),
                max_bonds_to_passivate,
                max_coordination_number - len(existing_vectors),
            )
            return best_missing[:num_bonds_to_add].tolist()

        return []


class BondDirectionForElement:
    bond_directions: BondDirections
    element: str


class BondDirectionForElementList(List[BondDirectionForElement]):
    values: List[BondDirectionForElement]
