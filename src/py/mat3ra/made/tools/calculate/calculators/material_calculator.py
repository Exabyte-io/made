from mat3ra.made.material import Material
from pydantic import BaseModel

from .material_calculator_parameters import MaterialCalculatorParameters


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
