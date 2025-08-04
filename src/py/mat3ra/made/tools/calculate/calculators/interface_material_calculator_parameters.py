from .material_calculator_parameters import MaterialCalculatorParameters


class InterfaceMaterialCalculatorParameters(MaterialCalculatorParameters):
    """
    Parameters specific to the calculation of interaction energies between
    an interface material's film and substrate.

    Args:
        shadowing_radius (float): Radius used to determine the surface atoms of the film or substrate
                                  for interaction calculations. Default is 2.5 Å.
    """

    shadowing_radius: float = 2.5
