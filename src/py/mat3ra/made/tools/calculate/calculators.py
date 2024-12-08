from typing import Callable

import numpy as np
from mat3ra.made.material import Material
from pydantic import BaseModel

from ..analyze.other import get_surface_atom_indices
from ..convert.utils import InterfacePartsEnum
from ..enums import SurfaceTypes
from ..modify import interface_get_part
from .interaction_functions import sum_of_inverse_distances_squared


class MaterialCalculatorParameters(BaseModel):
    """
    Defines the parameters for a material calculator.

    Args:
        interaction_function (Callable): A function used to calculate the interaction metric between
        sets of coordinates. The default function is sum_of_inverse_distances_squared.
    """

    interaction_function: Callable = sum_of_inverse_distances_squared


class InterfaceMaterialCalculatorParameters(MaterialCalculatorParameters):
    """
    Parameters specific to the calculation of interaction energies between
    an interface material's film and substrate.

    Args:
        shadowing_radius (float): Radius used to determine the surface atoms of the film or substrate
                                  for interaction calculations. Default is 2.5 Å.
    """

    shadowing_radius: float = 2.5


class MaterialCalculator(BaseModel):
    """
    A base class for performing calculations on materials.

    This class uses the parameters defined in MaterialCalculatorParameters to calculate
    interaction metrics between atoms or sets of coordinates within the material.

    Args:
        calculator_parameters (MaterialCalculatorParameters): Parameters controlling the calculator,
        including the interaction function.
    """

    calculator_parameters: MaterialCalculatorParameters = MaterialCalculatorParameters()

    def get_energy(self, material: Material):
        """
        Calculate the energy (or other metric) for a material.

        Args:
            material (Material): The material to calculate the interaction energy for.

        Returns:
            float: The interaction energy between the coordinates of the material,
            calculated using the specified interaction function.
        """
        return self.calculator_parameters.interaction_function(material.coordinates, material.coordinates)


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
