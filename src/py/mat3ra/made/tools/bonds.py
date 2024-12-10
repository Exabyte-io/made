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
        # Convert input_array to an ndarray of the same subclass
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

    def find_missing_directions(
        self,
        templates: List[np.ndarray],
        angle_tolerance: float = 0.1,
        max_bonds_to_add: int = 1,
    ) -> "BondDirections":
        """
        Reconstruct missing bonds for the atom based on the bond template best matching with the existing bonds.

        Args:
            templates (List[np.ndarray]): List of bond templates to match against.
            angle_tolerance (float): Tolerance for comparing angles between bond vectors.
            max_bonds_to_add (int): Maximum number of bonds to add.

        Returns:
            List[List[float]]: List of reconstructed bond vectors.
        """
        max_coordination_number = len(templates)

        if len(self) >= max_coordination_number:
            return BondDirections([])

        best_missing = None
        best_match_count = -1

        for template in templates:
            if self.size == 0:
                match_count = 0
            else:
                # TODO: optimize
                dot_matrix = np.dot(template, self.T)
                cosine_matrix = dot_matrix / (np.linalg.norm(template, axis=1)[:, None] * np.linalg.norm(self, axis=1))
                angles_matrix = np.arccos(np.clip(cosine_matrix, -1.0, 1.0))

                matches = np.any(angles_matrix < angle_tolerance, axis=1)
                match_count = np.sum(matches)

            missing = template[~matches] if self.size != 0 else template

            if match_count > best_match_count:
                best_match_count = match_count
                best_missing = missing

        if best_missing is not None:
            num_bonds_to_add = min(
                len(best_missing),
                max_bonds_to_add,
                max_coordination_number - len(self),
            )
            return BondDirections(best_missing[:num_bonds_to_add].tolist())

        return BondDirections([])


class BondDirectionsTemplatesForElement(BaseModel):
    bond_directions_templates: List[BondDirections]
    element: str

    model_config = ConfigDict(arbitrary_types_allowed=True)

    def to_ndarray(self):
        return np.array(self.bond_directions_templates)


class BondDirectionsForElementList(List[BondDirectionsTemplatesForElement]):
    pass
