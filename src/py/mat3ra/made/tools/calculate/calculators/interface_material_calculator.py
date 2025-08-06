from typing import Callable

import numpy as np
from mat3ra.made.material import Material
from mat3ra.made.tools.build.processed_structures.two_dimensional.passivation.enums import SurfaceTypesEnum

from ...analyze.other import get_surface_atom_indices
from ...convert.interface_parts_enum import InterfacePartsEnum
from ...modify import interface_get_part
from ..interaction_functions import sum_of_inverse_distances_squared
from .interface_material_calculator_parameters import InterfaceMaterialCalculatorParameters
from .material_calculator import MaterialCalculator


class InterfaceMaterialCalculator(MaterialCalculator):
    """
    A specialized calculator for computing the interaction energy between a film and a substrate
    in an interface material.

    This class extends MaterialCalculator and uses additional parameters specific to interface materials,
    such as the shadowing radius to detect surface atoms.

    Args:
        calculator_parameters (InterfaceMaterialCalculatorParameters): Parameters that include the
        shadowing radius and interaction function.
    """

    calculator_parameters: InterfaceMaterialCalculatorParameters = InterfaceMaterialCalculatorParameters()

    def get_energy(
        self,
        material: Material,
        shadowing_radius: float = 2.5,
        interaction_function: Callable = sum_of_inverse_distances_squared,
    ) -> float:
        """
        Calculate the interaction energy between the film and substrate in an interface material.

        This method uses the shadowing radius to detect surface atoms and applies the given
        interaction function to calculate the interaction between the film and substrate.

        Args:
            material (Material): The interface Material object consisting of both the film and substrate.
            shadowing_radius (float): The radius used to detect surface atoms for the interaction
                                      calculation. Defaults to 2.5 Å.
            interaction_function (Callable): A function to compute the interaction between the film and
                                             substrate. Defaults to sum_of_inverse_distances_squared.

        Returns:
            float: The calculated interaction energy between the film and substrate.
        """
        film_material = interface_get_part(material, part=InterfacePartsEnum.FILM)
        substrate_material = interface_get_part(material, part=InterfacePartsEnum.SUBSTRATE)

        film_surface_atom_indices = get_surface_atom_indices(
            film_material, SurfaceTypesEnum.BOTTOM, shadowing_radius=shadowing_radius
        )
        substrate_surface_atom_indices = get_surface_atom_indices(
            substrate_material, SurfaceTypesEnum.TOP, shadowing_radius=shadowing_radius
        )

        film_surface_atom_coordinates = film_material.basis.coordinates
        film_surface_atom_coordinates.filter_by_ids(film_surface_atom_indices)
        substrate_surface_atom_coordinates = substrate_material.basis.coordinates
        substrate_surface_atom_coordinates.filter_by_ids(substrate_surface_atom_indices)

        film_coordinates_values = np.array(film_surface_atom_coordinates.values)
        substrate_coordinates_values = np.array(substrate_surface_atom_coordinates.values)

        return interaction_function(film_coordinates_values, substrate_coordinates_values)
