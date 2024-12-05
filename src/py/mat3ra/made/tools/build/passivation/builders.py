from typing import Dict, List, Set
from mat3ra.made.material import Material
from pydantic import BaseModel, Field
import numpy as np
from ...enums import SurfaceTypes
from ...analyze import (
    get_surface_atom_indices,
    get_nearest_neighbors_vectors,
)
from ...modify import translate_to_z_level
from ...build import BaseBuilder
from .configuration import (
    PassivationConfiguration,
)


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
    Parameters for the CoordinationPassivationBuilder.

    Args:
        coordination_threshold (int): The coordination number threshold for atom to be considered undercoordinated.
        bonds_to_passivate (int): The maximum number of bonds to passivate for each undercoordinated atom.
        symmetry_tolerance (float): The tolerance for symmetry comparison of vectors for bonds.
    """

    coordination_threshold: int = Field(
        3, description="The coordination number threshold for an atom to be considered undercoordinated."
    )
    bonds_to_passivate: int = Field(
        1, description="The maximum number of bonds to passivate for each undercoordinated atom."
    )
    symmetry_tolerance: float = Field(0.1, description="The tolerance for symmetry comparison of vectors for bonds.")


class CoordinationBasedPassivationBuilder(PassivationBuilder):
    """
    Builder for passivating material surfaces based on coordination number analysis.

        This builder analyzes atomic coordination environments and reconstructs missing bonds
        using templates derived from fully coordinated atoms. It works by:

        1. Finding undercoordinated atoms that need passivation
        2. Creating bond vector templates for each element type by:
            - Collecting vectors from atoms with the highest valid coordination
            - Grouping similar vector arrangements into unique templates
        3. Reconstructing missing bonds by:
            - Matching existing bonds against templates
            - Finding the best-matching template for each atom
            - Adding missing bond vectors from the template
        4. Placing passivant atoms (typically H) along the reconstructed bond vectors
    """

    build_parameters: CoordinationBasedPassivationBuilderParameters
    _DefaultBuildParameters = CoordinationBasedPassivationBuilderParameters()
    _ConfigurationType = PassivationConfiguration

    def create_passivated_material(self, configuration: PassivationConfiguration) -> Material:
        """
        Create a passivated material based on the configuration.
        """
        material = super().create_passivated_material(configuration)
        coordination_analyzer = CoordinationAnalyzer(
            cutoff=self.build_parameters.shadowing_radius,
            coordination_threshold=self.build_parameters.coordination_threshold,
        )
        undercoordinated_atoms_indices = coordination_analyzer.get_undercoordinated_atom_indices(material)
        nearest_neighbors_vectors_array = get_nearest_neighbors_vectors(
            material=material, cutoff=self.build_parameters.shadowing_radius
        )
        templates = coordination_analyzer.find_template_vectors(
            nearest_neighbors_vectors_array.values, material.basis.elements.values
        )
        reconstructed_bonds = coordination_analyzer.reconstruct_missing_bonds(
            nearest_neighbors_vectors_array.values, material.basis.elements.values, templates
        )

        passivant_coordinates_values = self._get_passivant_coordinates(
            material,
            configuration,
            undercoordinated_atoms_indices,
            nearest_neighbors_vectors_array.values,
            reconstructed_bonds,
        )
        return self._add_passivant_atoms(material, passivant_coordinates_values, configuration.passivant)

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

        Returns:
            list: Coordinates where passivants should be added.
        """
        passivant_coordinates = []

        for idx in undercoordinated_atoms_indices:
            if idx in reconstructed_bonds:
                for bond_vector in reconstructed_bonds[idx]:
                    bond_vector = np.array(bond_vector)
                    if np.linalg.norm(bond_vector) == 0:
                        continue  # Avoid division by zero
                    normalized_bond = bond_vector / np.linalg.norm(bond_vector) * configuration.bond_length
                    passivant_bond_vector_crystal = material.basis.cell.convert_point_to_crystal(normalized_bond)
                    passivant_coordinates.append(
                        material.basis.coordinates.get_element_value_by_index(idx) + passivant_bond_vector_crystal
                    )
            else:
                # Fallback to nearest neighbors average vector method if no reconstructed bonds found
                vectors = nearest_neighbors_coordinates[idx]
                if vectors:
                    bond_vector = -np.mean(vectors, axis=0)
                    if np.linalg.norm(bond_vector) == 0:
                        continue  # Avoid division by zero
                    bond_vector = bond_vector / np.linalg.norm(bond_vector) * configuration.bond_length
                    passivant_bond_vector_crystal = material.basis.cell.convert_point_to_crystal(bond_vector)
                    passivant_coordinates.append(
                        material.basis.coordinates.get_element_value_by_index(idx) + passivant_bond_vector_crystal
                    )

        return passivant_coordinates


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

    def get_unique_coordination_numbers(self, material: Material) -> Set[int]:
        """
        Get the unique coordination numbers for all atoms in the material.

        Args:
            material (Material): Material object.

        Returns:
            Set[int]: A set of unique coordination numbers present in the material.
        """
        coordination_numbers = self.get_coordination_numbers(material)
        return set(coordination_numbers.values())

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
