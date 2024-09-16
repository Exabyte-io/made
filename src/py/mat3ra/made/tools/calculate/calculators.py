from typing import Callable

import numpy as np
from mat3ra.made.material import Material
from pydantic import BaseModel

from ..analyze import get_surface_atom_indices
from ..convert.utils import InterfacePartsEnum
from ..enums import SurfaceTypes
from ..modify import get_interface_part
from .interaction_functions import sum_of_inverse_distances_squared


class MaterialCalculatorParameters(BaseModel):
    interaction_function: Callable = sum_of_inverse_distances_squared


class InterfaceMaterialCalculatorParameters(MaterialCalculatorParameters):
    shadowing_radius: float = 2.5


class MaterialCalculator(BaseModel):
    calculator_parameters: MaterialCalculatorParameters = MaterialCalculatorParameters()

    def get_energy(self, material: Material):
        return self.calculator_parameters.interaction_function(material.coordinates, material.coordinates)


class InterfaceMaterialCalculator(MaterialCalculator):
    def get_energy(
        self,
        material: Material,
        shadowing_radius: float = 2.5,
        interaction_function: Callable = sum_of_inverse_distances_squared,
    ) -> float:
        """
        Calculate the interaction metric between the film and substrate.
        Args:
            material (Material): The interface Material object.
            shadowing_radius (float): The shadowing radius to detect the surface atoms, in Angstroms.
            interaction_function (Callable): The metric function to use for the calculation of the interaction.

        Returns:
            float: The calculated norm.
        """
        film_material = get_interface_part(material, part=InterfacePartsEnum.FILM)
        substrate_material = get_interface_part(material, part=InterfacePartsEnum.SUBSTRATE)
        film_surface_atom_indices = get_surface_atom_indices(
            film_material, SurfaceTypes.BOTTOM, shadowing_radius=shadowing_radius
        )
        substrate_surface_atom_indices = get_surface_atom_indices(
            substrate_material, SurfaceTypes.TOP, shadowing_radius=shadowing_radius
        )

        film_surface_atom_coordinates = film_material.basis.coordinates
        film_surface_atom_coordinates.filter_by_ids(film_surface_atom_indices)
        substrate_surface_atom_coordinates = substrate_material.basis.coordinates
        substrate_surface_atom_coordinates.filter_by_ids(substrate_surface_atom_indices)

        film_coordinates_values = np.array(film_surface_atom_coordinates.values)
        substrate_coordinates_values = np.array(substrate_surface_atom_coordinates.values)

        return interaction_function(film_coordinates_values, substrate_coordinates_values)
