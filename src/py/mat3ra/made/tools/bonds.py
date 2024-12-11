from enum import Enum
from typing import List

import numpy as np
from pydantic import BaseModel, ConfigDict


class BondDirectionsTemplatesEnum(list, Enum):
    OCTAHEDRAL = [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]]
    TETRAHEDRAL = [[1, 1, 1], [-1, -1, 1], [1, -1, -1], [-1, 1, -1]]
    PLANAR = [[1, 0, 0], [-1, 2, 0], [-1, -2, 0]]
    LINEAR = [[1, 0, 0], [-1, 0, 0]]


class BondDirections(np.ndarray):
    def __new__(cls, input_array, *args, **kwargs):
        """
        Create a new instance of BondDirections from an input array.

        Args:
            input_array (array-like): Input data to initialize the array.

        Returns:
            BondDirections: A new BondDirections instance.
        """
        obj = np.asarray(input_array).view(cls)
        return obj

    def __eq__(self, other):
        if not isinstance(other, BondDirections):
            return False
        if self.shape != other.shape:
            return False  # Prevents shape mismatches
        return self.are_bond_directions_similar(self, other)

    @staticmethod
    def are_bond_directions_similar(directions1: np.ndarray, directions2: np.ndarray, tolerance: float = 0.1) -> bool:
        """
        Check if two bond direction templates are similar.

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

    def find_missing_directions(
        self,
        bond_templates: List[np.ndarray],
        angle_tolerance: float = 0.1,
        max_bonds_to_add: int = 1,
    ) -> "BondDirections":
        """
        Reconstruct missing bonds for the atom based on the bond template that matches most of the existing bonds.

        Args:
            bond_templates (List[np.ndarray]): List of bond templates to match against.
            angle_tolerance (float): Tolerance for comparing angles between bond vectors.
            max_bonds_to_add (int): Maximum number of bonds to add.

        Returns:
            BondDirections: List of reconstructed bond vectors.
        """
        max_coordination_number = max(len(template) for template in bond_templates)

        if len(self) >= max_coordination_number:
            return BondDirections([])

        best_missing_directions = None
        highest_match_count = -1

        bond_templates = self._flatten_bond_templates(bond_templates)

        for bond_template in bond_templates:
            # Check which existing bonds match the template bonds
            is_matched = [
                any(self._are_vectors_similar(existing_bond, candidate_bond, angle_tolerance) for existing_bond in self)
                for candidate_bond in bond_template
            ]

            current_match_count = sum(is_matched)
            missing_directions = [
                candidate_bond for candidate_bond, matched in zip(bond_template, is_matched) if not matched
            ]

            # Select the template with the highest match count
            if current_match_count > highest_match_count:
                highest_match_count = current_match_count
                best_missing_directions = missing_directions

        if best_missing_directions is not None:
            num_bonds_to_add = min(
                len(best_missing_directions),
                max_bonds_to_add,
                max_coordination_number - len(self),
            )
            return BondDirections(best_missing_directions[:num_bonds_to_add])

        return BondDirections([])

    @staticmethod
    def _are_vectors_similar(vector1: np.ndarray, vector2: np.ndarray, tolerance: float = 0.1) -> bool:
        """
        Check if two bond vectors are similar.

        Args:
            vector1 (np.ndarray): First bond vector.
            vector2 (np.ndarray): Second bond vector.
            tolerance (float): Angle tolerance for comparison.

        Returns:
            bool: True if the bond vectors are similar, False otherwise.
        """
        dot_product = np.dot(vector1, vector2)
        norms = np.linalg.norm(vector1) * np.linalg.norm(vector2)
        cosine_angle = dot_product / norms if norms != 0 else 0.0
        angle = np.arccos(np.clip(cosine_angle, -1.0, 1.0))
        return angle < tolerance

    @staticmethod
    def _flatten_bond_templates(bond_templates: List[np.ndarray]) -> List[np.ndarray]:
        """
        Ensure that bond_templates is a list of np.ndarray with no extra nesting.

        Args:
            bond_templates (list): List of bond templates to preprocess.

        Returns:
            list: Flattened list of bond templates (np.ndarray).
        """
        flattened_templates = []
        for template in bond_templates:
            if isinstance(template, list) and all(isinstance(sub_template, np.ndarray) for sub_template in template):
                # If the template is a list of np.ndarray, extend the flattened list
                flattened_templates.extend(template)
            else:
                # Otherwise, add the template directly
                flattened_templates.append(template)
        return flattened_templates


class BondDirectionsTemplatesForElement(BaseModel):
    bond_directions_templates: List[BondDirections]
    element: str

    model_config = ConfigDict(arbitrary_types_allowed=True)

    def to_ndarray(self):
        return [np.array(template) for template in self.bond_directions_templates]


class BondDirectionsForElementList(List[BondDirectionsTemplatesForElement]):
    pass
