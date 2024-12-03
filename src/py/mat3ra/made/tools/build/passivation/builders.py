from typing import List

import numpy as np
from mat3ra.made.material import Material
from pydantic import BaseModel

from mat3ra.made.utils import ArrayWithIds
from ...enums import SurfaceTypes
from ...analyze import (
    get_surface_atom_indices,
    get_undercoordinated_atom_indices,
    get_nearest_neighbors_vectors,
    get_coordination_numbers,
)
from ...modify import translate_to_z_level
from ...build import BaseBuilder
from .configuration import (
    PassivationConfiguration,
)

import numpy as np
from collections import defaultdict
from typing import Dict, List, Tuple, Set
from dataclasses import dataclass


class PassivationBuilder(BaseBuilder):
    """
    Base class for passivation builders.
    """

    _GeneratedItemType = Material
    _ConfigurationType = PassivationConfiguration

    def _generate(self, configuration: BaseBuilder._ConfigurationType) -> List[Material]:
        return [self.create_passivated_material(configuration)]

    def _update_material_name(
        self, material: BaseBuilder._GeneratedItemType, configuration: BaseBuilder._ConfigurationType
    ) -> BaseBuilder._GeneratedItemType:
        material = super()._update_material_name(material, configuration)
        material.name += f" {configuration.passivant}-passivated"
        return material

    def create_passivated_material(self, configuration: BaseBuilder._ConfigurationType) -> Material:
        material = translate_to_z_level(configuration.slab, "center")
        return material

    def _add_passivant_atoms(self, material: Material, coordinates: list, passivant: str) -> Material:
        """
        Add passivant atoms to the provided coordinates in the material.

        Args:
            material (Material): The material object to add passivant atoms to.
            coordinates (list): The coordinates to add passivant atoms to.
            passivant (str): The chemical symbol of the passivating atom (e.g., 'H').

        Returns:
            Material: The material object with passivation atoms added.
        """
        for coord in coordinates:
            material.add_atom(passivant, coord)
        return material


class SurfacePassivationBuilderParameters(BaseModel):
    """
    Parameters for the SurfacePassivationBuilder, defining how atoms near the surface are
    detected and passivated.

    Args:
        shadowing_radius (float): Radius around each surface atom to exclude underlying atoms from passivation.
        depth (float): Depth from the topmost (or bottommost) atom into the material to consider for passivation,
                       accounting for features like islands, adatoms, and terraces.
    """

    shadowing_radius: float = 2.5
    depth: float = 5.0


class SurfacePassivationBuilder(PassivationBuilder):
    """
    Builder for passivating a surface.

    Detects surface atoms looking along Z axis and passivates either the top or bottom surface or both.
    """

    build_parameters: SurfacePassivationBuilderParameters
    _DefaultBuildParameters = SurfacePassivationBuilderParameters()
    _ConfigurationType = PassivationConfiguration

    def create_passivated_material(self, configuration: PassivationConfiguration) -> Material:
        material = super().create_passivated_material(configuration)
        passivant_coordinates_values = []

        if configuration.surface in (SurfaceTypes.BOTTOM, SurfaceTypes.BOTH):
            passivant_coordinates_values.extend(
                self._get_passivant_coordinates(material, SurfaceTypes.BOTTOM, configuration)
            )
        if configuration.surface in (SurfaceTypes.TOP, SurfaceTypes.BOTH):
            passivant_coordinates_values.extend(
                self._get_passivant_coordinates(material, SurfaceTypes.TOP, configuration)
            )

        return self._add_passivant_atoms(material, passivant_coordinates_values, configuration.passivant)

    def _get_passivant_coordinates(
        self, material: Material, surface: SurfaceTypes, configuration: PassivationConfiguration
    ):
        """
        Calculate the coordinates for placing passivants based on the specified surface type.

        Args:
            material (Material): Material to passivate.
            surface (SurfaceTypes): Surface type (TOP or BOTTOM).
            configuration (SurfacePassivationConfiguration): Configuration for passivation.

        Returns:
            list: Coordinates where passivants should be added.
        """
        surface_atoms_indices = get_surface_atom_indices(
            material, surface, self.build_parameters.shadowing_radius, self.build_parameters.depth
        )
        surface_atoms_coordinates = [
            material.basis.coordinates.get_element_value_by_index(i) for i in surface_atoms_indices
        ]
        bond_vector = (
            [0, 0, configuration.bond_length] if surface == SurfaceTypes.TOP else [0, 0, -configuration.bond_length]
        )
        passivant_bond_vector_crystal = material.basis.cell.convert_point_to_crystal(bond_vector)
        return (np.array(surface_atoms_coordinates) + np.array(passivant_bond_vector_crystal)).tolist()


class CoordinationBasedPassivationBuilderParameters(SurfacePassivationBuilderParameters):
    """
    Parameters for the  CoordinationPassivationBuilder.
    Args:
        coordination_threshold (int): The coordination number threshold for atom to be considered undercoordinated.
    """

    coordination_threshold: int = 3


class CoordinationBasedPassivationBuilder(PassivationBuilder):
    """
    Builder for passivating material based on coordination number of each atom.

    Detects atoms with coordination number below a threshold and passivates them.
    """

    _BuildParametersType = CoordinationBasedPassivationBuilderParameters
    _DefaultBuildParameters = CoordinationBasedPassivationBuilderParameters()

    def create_passivated_material(self, configuration: PassivationConfiguration) -> Material:
        material = super().create_passivated_material(configuration)
        all_indices = material.basis.coordinates.ids
        undercoordinated_atoms_indices = get_undercoordinated_atom_indices(
            material=material,
            indices=all_indices,
            cutoff=self.build_parameters.shadowing_radius,
            coordination_threshold=self.build_parameters.coordination_threshold,
        )
        nearest_neighbors_vectors_array = get_nearest_neighbors_vectors(
            material=material, indices=all_indices, cutoff=self.build_parameters.shadowing_radius
        )

        atom_vectors = self.convert_to_analysis_format(nearest_neighbors_vectors_array)
        atom_elements = {idx: material.basis.elements.get_element_value_by_index(idx) for idx in all_indices}
        reconstructed_bonds = self.analyze_crystal_structure(
            atom_vectors, atom_elements, symmetry_tolerance=0.1, distance_tolerance=0.1
        )

        passivant_coordinates_values = self._get_passivant_coordinates(
            material,
            configuration,
            undercoordinated_atoms_indices,
            nearest_neighbors_vectors_array.values,
            reconstructed_bonds,
        )
        return self._add_passivant_atoms(material, passivant_coordinates_values, configuration.passivant)

    def convert_to_analysis_format(self, array_with_ids: ArrayWithIds) -> Dict[int, List[List[float]]]:
        """
        Converts ArrayWithIds format to the format needed for crystal structure analysis.

        Expected input format:
        - array_with_ids.values contains list of tuples/lists where:
            - Each item is (vector_list, element)
            - vector_list is a list of 3D vectors

        Returns:
        - atom_vectors: Dict[int, List[List[float]]] - mapping of atom ID to its vectors
        """
        atom_vectors = {}

        for atom_id, vectors in zip(array_with_ids.ids, array_with_ids.values):
            # Convert vectors to list of lists if they're not already
            if isinstance(vectors, np.ndarray):
                vectors = vectors.tolist()
            elif not isinstance(vectors, list):
                vectors = [vectors]

            # Ensure vectors are in the correct format
            formatted_vectors = []
            for vec in vectors:
                if isinstance(vec, (list, tuple, np.ndarray)):
                    formatted_vectors.append([float(v) for v in vec])
                else:
                    raise ValueError(f"Invalid vector format for atom {atom_id}: {vec}")

            atom_vectors[atom_id] = formatted_vectors

        return atom_vectors

    @dataclass
    class AtomData:
        id: int
        element: str
        vectors: np.ndarray  # Shape: (N, 3) where N is number of vectors

    def align_vectors(self, vectors1: np.ndarray, vectors2: np.ndarray) -> Tuple[float, np.ndarray]:
        """
        Align two sets of vectors using Kabsch algorithm.
        Returns RMSD and rotation matrix.
        """
        # Center both vector sets
        centroid1 = np.mean(vectors1, axis=0)
        centroid2 = np.mean(vectors2, axis=0)
        vectors1_centered = vectors1 - centroid1
        vectors2_centered = vectors2 - centroid2

        # Calculate covariance matrix
        covariance_matrix = vectors1_centered.T @ vectors2_centered

        # SVD decomposition
        U, _, Vt = np.linalg.svd(covariance_matrix)

        # Calculate rotation matrix
        rotation_matrix = Vt.T @ U.T

        # Apply rotation
        aligned_vectors = vectors2_centered @ rotation_matrix.T + centroid1

        # Calculate RMSD
        rmsd = np.sqrt(np.mean(np.sum((vectors1 - aligned_vectors) ** 2, axis=1)))

        return rmsd, rotation_matrix

    def find_coordination_templates(
        self, atoms_data: List[AtomData], symmetry_tolerance: float = 0.1
    ) -> Dict[str, np.ndarray]:
        """
        Find ideal coordination templates for each element.
        Returns a dictionary of element -> template vectors.
        """
        # Group atoms by element
        element_groups = defaultdict(list)
        for atom in atoms_data:
            element_groups[atom.element].append(atom)

        templates = {}

        for element, atoms in element_groups.items():
            # Find atom with maximum number of vectors
            max_vectors_atom = max(atoms, key=lambda x: len(x.vectors))
            max_coords = len(max_vectors_atom.vectors)

            # Find all atoms with same number of vectors
            template_candidates = [atom for atom in atoms if len(atom.vectors) == max_coords]

            # If only one candidate, use its vectors as template
            if len(template_candidates) == 1:
                templates[element] = max_vectors_atom.vectors
                continue

            # Align all candidates and check symmetry
            reference_vectors = template_candidates[0].vectors
            aligned_sets = [reference_vectors]

            for candidate in template_candidates[1:]:
                rmsd, rotation = self.align_vectors(reference_vectors, candidate.vectors)
                if rmsd <= symmetry_tolerance:
                    aligned_vectors = candidate.vectors @ rotation.T
                    aligned_sets.append(aligned_vectors)

            # Average the aligned vector sets
            templates[element] = np.mean(aligned_sets, axis=0)

        return templates

    def reconstruct_missing_bonds(
        self, atoms_data: List[AtomData], templates: Dict[str, np.ndarray], distance_tolerance: float = 0.1
    ) -> Dict[int, np.ndarray]:
        """
        Reconstruct missing bonds based on templates, accounting for existing bonds.
        Returns dictionary of atom_id -> reconstructed vectors.
        """
        missing_bonds = {}

        for atom in atoms_data:
            template = templates[atom.element]
            expected_coords = len(template)
            current_coords = len(atom.vectors)

            # Only proceed if we have fewer bonds than expected
            if current_coords < expected_coords:
                # For atoms with existing bonds, align template with them
                if current_coords > 0:
                    # Use existing bonds to align template
                    rmsd, rotation = self.align_vectors(template[:current_coords], atom.vectors)
                    template_rotated = template @ rotation.T
                else:
                    # No existing bonds, use template as is
                    template_rotated = template

                # Check each template vector against existing vectors
                existing_vectors = atom.vectors
                missing = []

                for template_vector in template_rotated:
                    # Flag to check if this template vector matches any existing vector
                    is_new_direction = True

                    for existing_vector in existing_vectors:
                        # Calculate angular difference between vectors
                        cos_angle = np.dot(template_vector, existing_vector) / (
                            np.linalg.norm(template_vector) * np.linalg.norm(existing_vector)
                        )
                        cos_angle = np.clip(cos_angle, -1.0, 1.0)  # Handle numerical errors
                        angle = np.arccos(cos_angle)

                        # If vectors are nearly parallel or antiparallel, consider them matching
                        if angle < distance_tolerance or abs(angle - np.pi) < distance_tolerance:
                            is_new_direction = False
                            break

                    if is_new_direction:
                        missing.append(template_vector)

                if missing:
                    missing_bonds[atom.id] = np.array(missing)

        return missing_bonds

    def analyze_crystal_structure(
        self,
        atom_vectors: Dict[int, List[List[float]]],
        atom_elements: Dict[int, str],
        symmetry_tolerance: float = 0.1,
        distance_tolerance: float = 0.1,
    ) -> Dict[int, List[List[float]]]:
        """
        Main function to analyze crystal structure and reconstruct missing bonds.

        Parameters:
        - atom_vectors: Dictionary of atom_id -> list of vectors
        - atom_elements: Dictionary of atom_id -> chemical element
        - symmetry_tolerance: Maximum RMSD for considering vectors symmetric
        - distance_tolerance: Maximum distance for matching vectors

        Returns:
        - Dictionary of atom_id -> list of reconstructed vectors
        """
        # Convert input data to internal format
        atoms_data = [
            self.AtomData(id=atom_id, element=atom_elements[atom_id], vectors=np.array(vectors))
            for atom_id, vectors in atom_vectors.items()
        ]

        # Find coordination templates
        templates = self.find_coordination_templates(atoms_data, symmetry_tolerance=symmetry_tolerance)

        # Reconstruct missing bonds
        missing_bonds = self.reconstruct_missing_bonds(atoms_data, templates, distance_tolerance=distance_tolerance)

        # Convert numpy arrays back to lists for return
        return {atom_id: vectors.tolist() for atom_id, vectors in missing_bonds.items()}

    def _get_passivant_coordinates(
        self,
        material: Material,
        configuration: PassivationConfiguration,
        undercoordinated_atoms_indices: list,
        nearest_neighbors_coordinates: list,
        reconstructed_bonds: Dict[int, List[List[float]]],
    ):
        """
        Calculate the coordinates for placing passivating atoms based on reconstructed bonds.

        Args:
            material (Material): Material to passivate.
            configuration (PassivationConfiguration): Configuration for passivation.
            undercoordinated_atoms_indices (list): Indices of undercoordinated atoms.
            nearest_neighbors_coordinates (list): List of nearest neighbor vectors for each atom.
            reconstructed_bonds (Dict[int, List[List[float]]]): Dictionary of reconstructed bonds for each atom.
        """
        passivant_coordinates = []

        for idx in undercoordinated_atoms_indices:
            if idx in reconstructed_bonds:
                # Use reconstructed bonds for positioning passivants
                for bond_vector in reconstructed_bonds[idx]:
                    # Normalize bond vector and scale to desired bond length
                    bond_vector = np.array(bond_vector)
                    normalized_bond = bond_vector / np.linalg.norm(bond_vector) * configuration.bond_length

                    # Convert to crystal coordinates
                    passivant_bond_vector_crystal = material.basis.cell.convert_point_to_crystal(normalized_bond)

                    # Add passivant coordinate
                    passivant_coordinates.append(
                        material.basis.coordinates.get_element_value_by_index(idx) + passivant_bond_vector_crystal
                    )
            else:
                # Fallback to original method if no reconstructed bonds (shouldn't happen for undercoordinated atoms)
                vectors = nearest_neighbors_coordinates[idx]
                bond_vector = -np.mean(vectors, axis=0)
                bond_vector = bond_vector / np.linalg.norm(bond_vector) * configuration.bond_length
                passivant_bond_vector_crystal = material.basis.cell.convert_point_to_crystal(bond_vector)

                passivant_coordinates.append(
                    material.basis.coordinates.get_element_value_by_index(idx) + passivant_bond_vector_crystal
                )

        return passivant_coordinates

    def get_unique_coordination_numbers(self, material: Material) -> set:
        """
        Get unique coordination numbers for all atoms in the material for current builder parameters.

        Args:
            material (Material): The material object.

        Returns:
            set: The coordination numbers for all atoms in the material.
        """

        coordination_numbers = set(
            get_coordination_numbers(material=material, cutoff=self.build_parameters.shadowing_radius).values
        )
        return coordination_numbers
