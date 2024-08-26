from typing import List
from ...build import BaseBuilder
from mat3ra.made.material import Material
from .configuration import PassivationConfiguration


class PassivationBuilder(BaseBuilder):
    """
    Base class for passivation builders.
    """

    def _generate(self, configuration: PassivationConfiguration) -> List[Material]:
        return [self.create_passivated_material(configuration)]

    def _update_material_name(self, material, configuration):
        material = super()._update_material_name(material, configuration)
        material.name += f" {configuration.passivant}-passivated"
        return material

    def _update_material_basis(self, material, configuration):
        return super()._update_material_basis(material, configuration)

    def create_passivated_material(self, configuration: PassivationConfiguration) -> Material:
        raise NotImplementedError

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


class SurfacePassivationBuilder(PassivationBuilder):
    """
    Builder for passivating a surface.

    Detects surface atoms lokking along Z axis and passivates either the top or bottom surface or both.
    """

    def create_passivated_material(self, configuration: PassivationConfiguration) -> Material:
        # Reference to the passivate_surface function
        # startLine: 56
        # endLine: 99
        pass


class EdgePassivationBuilder(PassivationBuilder):
    """
    Builder for passivating an edge.

    Detects edge atoms looking perpendicular to the Z axis and passivates them.
    """

    def create_passivated_material(self, configuration: PassivationConfiguration) -> Material:
        # Reference to the passivate_edges function
        # startLine: 100
        # endLine: 143
        pass
