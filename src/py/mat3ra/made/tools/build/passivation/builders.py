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

        atom_vectors = nearest_neighbors_vectors_array
        atom_elements = material.basis.elements
        symmetry_tolerance = 0.1

        reconstructed_bonds = self.reconstruct_missing_bonds(atom_vectors, atom_elements, symmetry_tolerance)
        # Convert numpy arrays back to lists for return

        passivant_coordinates_values = self._get_passivant_coordinates(
            material,
            configuration,
            undercoordinated_atoms_indices,
            nearest_neighbors_vectors_array.values,
            reconstructed_bonds,
        )
        return self._add_passivant_atoms(material, passivant_coordinates_values, configuration.passivant)

    def vectors_are_similar(self, vec1: np.ndarray, vec2: np.ndarray, tolerance: float = 0.1) -> bool:
        """Check if two vectors are similar within tolerance."""
        v1_norm = vec1 / np.linalg.norm(vec1)
        v2_norm = vec2 / np.linalg.norm(vec2)
        angle = np.arccos(np.clip(np.dot(v1_norm, v2_norm), -1.0, 1.0))
        return angle < tolerance

    def templates_are_similar(self, template1: np.ndarray, template2: np.ndarray, tolerance: float = 0.1) -> bool:
        """Check if two templates are similar (have similar vectors)."""
        if len(template1) != len(template2):
            return False

        # Check if all vectors in template1 have a match in template2
        unmatched = list(range(len(template2)))
        for v1 in template1:
            found_match = False
            for i in unmatched:
                if self.vectors_are_similar(v1, template2[i], tolerance):
                    unmatched.remove(i)
                    found_match = True
                    break
            if not found_match:
                return False
        return True

    def find_template_vectors(
        self, atom_vectors: ArrayWithIds, atom_elements: ArrayWithIds
    ) -> Dict[str, List[np.ndarray]]:
        """Find unique template vector sets with maximum coordination for each element."""
        element_templates = {}

        for element in set(atom_elements.values):
            # Get all vectors for this element
            element_indices = [i for i, e in enumerate(atom_elements.values) if e == element]
            element_vector_lists = [np.array(atom_vectors.values[i]) for i in element_indices]

            # Find max coordination
            max_coord = max(len(vectors) for vectors in element_vector_lists)
            max_coord_vectors = [v for v in element_vector_lists if len(v) == max_coord]

            # Group similar templates
            unique_templates: List[np.ndarray] = []
            for template in max_coord_vectors:
                # Check if similar template exists
                template_exists = False
                for existing in unique_templates:
                    if self.templates_are_similar(template, existing):
                        template_exists = True
                        break
                if not template_exists:
                    unique_templates.append(template)

            element_templates[element] = unique_templates

        return element_templates

    def reconstruct_missing_bonds(
        self, atom_vectors: ArrayWithIds, atom_elements: ArrayWithIds, angle_tolerance: float = 0.1
    ) -> Dict[int, List[List[float]]]:
        """Find missing bonds using template matching."""
        templates = self.find_template_vectors(atom_vectors, atom_elements)
        missing_bonds = {}

        for idx, (vectors, element) in enumerate(zip(atom_vectors.values, atom_elements.values)):
            if element not in templates:
                continue

            # Skip if we have max coordination
            if len(vectors) >= len(templates[element][0]):
                continue

            existing_vectors = np.array(vectors) if vectors else np.empty((0, 3))

            # Try each template to find best match
            best_missing = None
            best_match_count = -1

            for template in templates[element]:
                current_missing = []
                current_match_count = 0

                for template_v in template:
                    # Check if template vector matches any existing vector
                    has_match = False
                    for existing_v in existing_vectors:
                        if self.vectors_are_similar(template_v, existing_v, angle_tolerance):
                            has_match = True
                            current_match_count += 1
                            break

                    if not has_match:
                        current_missing.append(template_v)

                if current_match_count > best_match_count:
                    best_match_count = current_match_count
                    best_missing = current_missing

            if best_missing:
                missing_bonds[atom_vectors.ids[idx]] = best_missing

        return {k: [v.tolist() for v in vectors] for k, vectors in missing_bonds.items()}

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
                for bond_vector in reconstructed_bonds[idx]:
                    bond_vector = np.array(bond_vector)
                    normalized_bond = bond_vector / np.linalg.norm(bond_vector) * configuration.bond_length

                    passivant_bond_vector_crystal = material.basis.cell.convert_point_to_crystal(normalized_bond)

                    passivant_coordinates.append(
                        material.basis.coordinates.get_element_value_by_index(idx) + passivant_bond_vector_crystal
                    )
            # else:
            #     # Fallback to original method if no reconstructed bonds (shouldn't happen for undercoordinated atoms)
            #     vectors = nearest_neighbors_coordinates[idx]
            #     bond_vector = -np.mean(vectors, axis=0)
            #     bond_vector = bond_vector / np.linalg.norm(bond_vector) * configuration.bond_length
            #     passivant_bond_vector_crystal = material.basis.cell.convert_point_to_crystal(bond_vector)
            #
            #     passivant_coordinates.append(
            #         material.basis.coordinates.get_element_value_by_index(idx) + passivant_bond_vector_crystal
            #     )

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
