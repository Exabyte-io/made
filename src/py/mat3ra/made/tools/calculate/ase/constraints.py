from typing import List, Optional

import numpy as np
from pydantic import BaseModel

from ...third_party import ASEAtoms


class InterfaceConstraint(BaseModel):
    """
    Based on  https://wiki.fysik.dtu.dk/ase/ase/constraints.html#making-your-own-constraint-class
    """

    class Config:
        arbitrary_types_allowed = True

    def adjust_positions(self, atoms: ASEAtoms, new_positions: np.ndarray) -> None:
        """
        Adjust the positions of the atoms in the system.

        Args:
            atoms (ASEAtoms): The ASE Atoms object.
            new_positions (numpy.ndarray): The new positions to adjust in place.
        """
        raise NotImplementedError

    def adjust_forces(self, forces: np.ndarray) -> None:
        """
        Adjust the forces on the atoms in the system.

        Args:
            forces (numpy.ndarray): The forces array to adjust in place.
        """
        raise NotImplementedError


class RigidFilmXYInterfaceConstraint(InterfaceConstraint):
    """
    Custom constraint to allow only rigid translation in x and y for film atoms.
    """

    film_indices: List[int] = []
    initial_positions: Optional[np.ndarray] = None
    reference_center: Optional[np.ndarray] = None

    def adjust_positions(self, atoms: ASEAtoms, new_positions: np.ndarray) -> None:
        """
        Adjust the positions of film atoms to maintain rigidity in x and y.

        Args:
            atoms (ase.Atoms): The ASE Atoms object.
            new_positions (numpy.ndarray): The new positions to adjust in place.
        """
        if self.initial_positions is None:
            self.initial_positions = atoms.positions[self.film_indices].copy()
            self.reference_center = self.initial_positions.mean(axis=0)

        current_center = new_positions[self.film_indices].mean(axis=0)
        desired_translation = current_center - self.reference_center
        desired_translation[2] = 0.0

        # Update positions to maintain rigid body translation in x and y
        for i, idx in enumerate(self.film_indices):
            new_positions[idx, 0] = self.initial_positions[i, 0] + desired_translation[0]
            new_positions[idx, 1] = self.initial_positions[i, 1] + desired_translation[1]

    def adjust_forces(self, forces: np.ndarray) -> None:
        """
        Adjust the forces on film atoms to ensure rigidity.

        Args:
            forces (numpy.ndarray): The forces array to adjust in place.
        """
        if self.initial_positions is None:
            return  # No adjustment needed if initial positions are not set

        # Project forces onto x and y for collective translation
        net_force_x = np.sum(forces[self.film_indices, 0])
        net_force_y = np.sum(forces[self.film_indices, 1])

        # Distribute the net force equally to all film atoms to maintain rigidity
        num_film_atoms = len(self.film_indices)
        if num_film_atoms == 0:
            return  # Avoid division by zero

        distributed_force_x = net_force_x / num_film_atoms
        distributed_force_y = net_force_y / num_film_atoms

        for idx in self.film_indices:
            forces[idx, 0] = distributed_force_x
            forces[idx, 1] = distributed_force_y
