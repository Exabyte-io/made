from typing import List
from mat3ra.made.material import Material
from mat3ra.made.tools.bonds import BondDirections
from pydantic import BaseModel
import numpy as np

from ...analyze.material import MaterialWithCrystalSites
from ...enums import SurfaceTypes
from ...analyze.other import get_surface_atom_indices
from ...modify import translate_to_z_level
from ...build import BaseBuilder
from .configuration import PassivationConfiguration


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
        material = configuration.slab.clone()
        # material.to_cartesian()
        material = translate_to_z_level(material, "center")
        # material.to_crystal()
        return material

    def _add_passivant_atoms(
        self, material: Material, coordinates: list, passivant: str, use_cartesian_coordinates=False
    ) -> Material:
        """
        Add passivant atoms to the provided coordinates in the material.

        Args:
            material (Material): The material object to add passivant atoms to.
            coordinates (list): The coordinates to add passivant atoms to.
            passivant (str): The chemical symbol of the passivating atom (e.g., 'H').
            use_cartesian_coordinates (bool): Whether the coordinates are in Cartesian units (or crystal by default).

        Returns:
            Material: The material object with passivation atoms added.
        """
        for coord in coordinates:
            material.add_atom(passivant, coord, use_cartesian_coordinates)
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

    coordination_threshold: int = 3  # The coordination number threshold for an atom to be considered undercoordinated.
    bonds_to_passivate: int = 1  # The maximum number of bonds to passivate for each undercoordinated atom.
    symmetry_tolerance: float = 0.1  # The tolerance for symmetry comparison of vectors for bonds.


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
        material_with_crystal_sites = MaterialWithCrystalSites.from_material(material)
        material_with_crystal_sites.analyze()

        undercoordinated_atoms_indices = material_with_crystal_sites.get_undercoordinated_atom_indices(
            cutoff=self.build_parameters.shadowing_radius,
            coordination_threshold=self.build_parameters.coordination_threshold,
        )
        # TODO: bonds_templates will be passed from the configuration in the "controlled" version of this class
        bonds_templates = material_with_crystal_sites.find_unique_bond_directions()
        reconstructed_bonds = material_with_crystal_sites.find_missing_bonds_for_all_sites(
            bond_templates_list=bonds_templates,
            max_bonds_to_add=self.build_parameters.bonds_to_passivate,
            angle_tolerance=self.build_parameters.symmetry_tolerance,
        )
        passivant_coordinates_values = self._get_passivant_coordinates(
            material,
            configuration,
            undercoordinated_atoms_indices,
            reconstructed_bonds,
        )
        return self._add_passivant_atoms(material, passivant_coordinates_values, configuration.passivant, False)

    def _get_passivant_coordinates(
        self,
        material: Material,
        configuration: PassivationConfiguration,
        undercoordinated_atoms_indices: list,
        reconstructed_bonds: List[BondDirections],
    ):
        """
        Calculate the coordinates for placing passivating atoms based on reconstructed bonds.

        Args:
            material (Material): Material to passivate.
            configuration (PassivationConfiguration): Configuration for passivation.
            undercoordinated_atoms_indices (list): Indices of undercoordinated atoms.
            reconstructed_bonds (List[BondDirections]): Reconstructed bond directions.

        Returns:
            list: Coordinates where passivants should be added.
        """
        passivant_coordinates = []

        for idx in undercoordinated_atoms_indices:
            for bond_vector in reconstructed_bonds[idx]:
                bond_vector_np = np.array(bond_vector)
                if np.linalg.norm(bond_vector_np) == 0:
                    continue  # Avoid division by zero
                normalized_bond = bond_vector_np / np.linalg.norm(bond_vector_np) * configuration.bond_length
                normalized_bond_crystal = material.basis.cell.convert_point_to_crystal(normalized_bond)
                passivant_coordinates.append(
                    np.array(material.basis.coordinates.get_element_value_by_index(idx))
                    + np.array(normalized_bond_crystal)
                )

        return passivant_coordinates
