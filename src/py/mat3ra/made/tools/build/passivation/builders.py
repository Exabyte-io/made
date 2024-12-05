from typing import Dict, List
from mat3ra.made.material import Material
from pydantic import BaseModel
import numpy as np
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
        bonds_to_passivate (int): The maximum number of bonds to passivate for each undercoordinated atom.
        symmetry_tolerance (float): The tolerance for symmetry comparison of vectors for bonds.
    """

    coordination_threshold: int = 3
    bonds_to_passivate: int = 1
    symmetry_tolerance: float = 0.1


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

    _BuildParametersType = CoordinationBasedPassivationBuilderParameters
    _DefaultBuildParameters = CoordinationBasedPassivationBuilderParameters()

    def create_passivated_material(self, configuration: PassivationConfiguration) -> Material:
        """
        Create a passivated material based on the configuration.
        """
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

        reconstructed_bonds = self.reconstruct_missing_bonds(
            nearest_neighbor_vectors_array=nearest_neighbors_vectors_array,
            chemical_elements_array=material.basis.elements,
            max_bonds_to_passivate=self.build_parameters.bonds_to_passivate,
            angle_tolerance=self.build_parameters.symmetry_tolerance,
        )

        passivant_coordinates_values = self._get_passivant_coordinates(
            material,
            configuration,
            undercoordinated_atoms_indices,
            nearest_neighbors_vectors_array.values,
            reconstructed_bonds,
        )
        return self._add_passivant_atoms(material, passivant_coordinates_values, configuration.passivant)

    @staticmethod
    def vectors_are_similar(vec1: np.ndarray, vec2: np.ndarray, tolerance: float = 0.1) -> bool:
        """Check if two vectors are similar within tolerance."""
        dot_product = np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))
        angle = np.arccos(np.clip(dot_product, -1.0, 1.0))
        return angle < tolerance

    @staticmethod
    def are_bonds_templates_similar(template1: np.ndarray, template2: np.ndarray, tolerance: float = 0.1) -> bool:
        """Check if two templates are similar."""
        if len(template1) != len(template2):
            return False

        # Calculate pairwise angles between vectors in both templates
        dot_matrix = np.dot(template1, template2.T)
        norms1 = np.linalg.norm(template1, axis=1)
        norms2 = np.linalg.norm(template2, axis=1)
        cosine_matrix = dot_matrix / np.outer(norms1, norms2)
        angles_matrix = np.arccos(np.clip(cosine_matrix, -1.0, 1.0))

        # Check if every vector in template1 matches a vector in template2
        unmatched = list(range(len(template2)))
        for i, angle_row in enumerate(angles_matrix):
            matches = np.where(angle_row < tolerance)[0]
            if len(matches) == 0:
                return False
            unmatched.remove(matches.tolist()[0])

        return True

    def find_template_vectors(
        self, atom_vectors: ArrayWithIds, atom_elements: ArrayWithIds
    ) -> Dict[str, List[np.ndarray]]:
        """
        Find unique templates (set of vectors to the nearest neighbor for fully coordinated atoms) for each element type

        Args:
            atom_vectors (ArrayWithIds): Array with atom vectors.
            atom_elements (ArrayWithIds): Array with atom elements.

        Returns:
            Dict[str, List[np.ndarray]]: Dictionary with element type as key and list of unique templates as value.
        """
        element_templates = {}

        for element in set(atom_elements.values):
            element_indices = [i for i, e in enumerate(atom_elements.values) if e == element]
            element_vector_lists = [np.array(atom_vectors.values[i]) for i in element_indices]

            max_coord = max(len(vectors) for vectors in element_vector_lists)
            max_coord_vectors = [v for v in element_vector_lists if len(v) == max_coord]

            unique_templates: List[np.ndarray] = []
            for template in max_coord_vectors:
                if not any(self.are_bonds_templates_similar(template, existing) for existing in unique_templates):
                    unique_templates.append(template)

            element_templates[element] = unique_templates

        return element_templates

    def reconstruct_missing_bonds(
        self,
        nearest_neighbor_vectors_array: ArrayWithIds,
        chemical_elements_array: ArrayWithIds,
        max_bonds_to_passivate: int = 2,
        angle_tolerance: float = 0.1,
    ) -> Dict[int, List[List[float]]]:
        """
        Reconstruct missing bonds for undercoordinated atoms based on templates and existing bonds.

        Args:
            nearest_neighbor_vectors_array (ArrayWithIds): Array with nearest neighbor vectors for each atom.
            chemical_elements_array (ArrayWithIds): Array with chemical elements for each atom.
            max_bonds_to_passivate (int): Maximum number of bonds to passivate for each undercoordinated atom.
            angle_tolerance (float): Tolerance for symmetry comparison of vectors for bonds.

        Returns:
            Dict[int, List[List[float]]]: Dict with atom indices as keys and list of missing bond vectors as values.
        """
        templates = self.find_template_vectors(nearest_neighbor_vectors_array, chemical_elements_array)
        missing_bonds = {}

        for idx, (vectors, element) in enumerate(
            zip(nearest_neighbor_vectors_array.values, chemical_elements_array.values)
        ):
            if element not in templates:
                continue
            existing_vectors = np.array(vectors) if vectors else np.empty((0, 3))
            max_coordination = len(templates[element][0])

            if len(existing_vectors) >= max_coordination:
                continue

            best_missing = None
            best_match_count = -1

            for template in templates[element]:
                dot_matrix = np.dot(template, existing_vectors.T)
                cosine_matrix = dot_matrix / (
                    np.linalg.norm(template, axis=1)[:, None] * np.linalg.norm(existing_vectors, axis=1)
                )
                angles_matrix = np.arccos(np.clip(cosine_matrix, -1.0, 1.0))

                matches = np.any(angles_matrix < angle_tolerance, axis=1)
                match_count = np.sum(matches)
                missing = template[~matches]

                if match_count > best_match_count:
                    best_match_count = match_count
                    best_missing = missing

            if best_missing is not None:
                num_bonds_to_add = min(
                    len(best_missing),
                    max_bonds_to_passivate,
                    max_coordination - len(existing_vectors),
                )
                missing_bonds[nearest_neighbor_vectors_array.ids[idx]] = best_missing[:num_bonds_to_add].tolist()

        return missing_bonds

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
                    normalized_bond = bond_vector / np.linalg.norm(bond_vector) * configuration.bond_length
                    passivant_bond_vector_crystal = material.basis.cell.convert_point_to_crystal(normalized_bond)
                    passivant_coordinates.append(
                        material.basis.coordinates.get_element_value_by_index(idx) + passivant_bond_vector_crystal
                    )
            else:
                # Fallback to nearest neighbors average vector method if no reconstructed bonds found
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
