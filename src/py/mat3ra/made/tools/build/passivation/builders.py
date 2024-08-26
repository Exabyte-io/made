from typing import List

import numpy as np
from mat3ra.made.material import Material
from pydantic import BaseModel

from .enums import SurfaceTypes
from ...analyze import get_surface_atoms_indices
from ...modify import translate_to_z_level
from ...build import BaseBuilder
from .configuration import PassivationConfiguration, SurfacePassivationConfiguration, EdgePassivationConfiguration


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
    Parameters for the SurfacePassivationBuilder.

    Args:
        shadowing_radius (float): Radius for atoms shadowing underlying from passivation, in Angstroms.
        depth (float): Depth from the top to look for exposed surface atoms to passivate, in Angstroms.
    """

    shadowing_radius: float = 2.5
    depth: float = 5.0


class SurfacePassivationBuilder(PassivationBuilder):
    """
    Builder for passivating a surface.

    Detects surface atoms looking along Z axis and passivates either the top or bottom surface or both.
    """

    def create_passivated_material(self, configuration: SurfacePassivationConfiguration) -> Material:
        material = super().create_passivated_material(configuration)
        passivant_coordinates_values_top = np.array([])
        passivant_coordinates_values_bottom = np.array([])

        if configuration.surface == SurfaceTypes.TOP or configuration.surface == SurfaceTypes.BOTH:
            passivant_coordinates_values_top = self._get_passivant_coordinates(
                material,
                SurfaceTypes.TOP,
                configuration.bond_length,
                self.build_parameters.shadowing_radius,
                self.build_parameters.depth,
            )

        if configuration.surface == SurfaceTypes.BOTTOM or configuration.surface == SurfaceTypes.BOTH:
            passivant_coordinates_values_bottom = self._get_passivant_coordinates(
                material,
                SurfaceTypes.BOTTOM,
                configuration.bond_length,
                self.build_parameters.shadowing_radius,
                self.build_parameters.depth,
            )

        passivant_coordinates_values = (
            passivant_coordinates_values_bottom.tolist() + passivant_coordinates_values_top.tolist()
        )
        return self._add_passivant_atoms(material, passivant_coordinates_values, configuration.passivant)

    def _get_passivant_coordinates(self, material, surface, bond_length, shadowing_radius, depth):
        surface_atoms_indices = get_surface_atoms_indices(
            material=material,
            surface=surface,
            shadowing_radius=shadowing_radius,
            depth=depth,
        )
        surface_atoms_coordinates = [
            material.basis.coordinates.get_element_value_by_index(i) for i in surface_atoms_indices
        ]
        bond_vector = [0, 0, bond_length] if surface == SurfaceTypes.TOP else [0, 0, -bond_length]
        passivant_bond_vector_crystal = material.basis.cell.convert_point_to_crystal(bond_vector)
        return np.array(surface_atoms_coordinates) + np.array(passivant_bond_vector_crystal)


class EdgePassivationBuilder(PassivationBuilder):
    """
    Builder for passivating an edge.

    Detects edge atoms looking perpendicular to the Z axis and passivates them.
    """

    def create_passivated_material(self, configuration: EdgePassivationConfiguration) -> Material:
        pass
