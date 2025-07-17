import sympy as sp

from mat3ra.made.material import Material
from .configuration import (
    PerturbationConfiguration,
)
from .builders import PerturbationBuilder
from ...utils.perturbation import PerturbationFunctionHolder


def create_perturbation(
    material: Material,
    perturbation_function: sp.Expr,
    use_cartesian_coordinates: bool = True,
    is_isometric: bool = False,
):
    """
    Create a perturbation for a material.
    With a function of delta for coordinate [∆x, ∆y, ∆z] = f(x, y, z) applied to each atom in the material.

    Args:
        material (Material): The material to be perturbed.
        perturbation_function (PerturbationFunctionHolder): The perturbation function to apply.
        is_isometric (bool): If True, the lattice will be adjusted to preserve distances between atoms.

    Returns:
        Material: The perturbed material.
    """
    variables = [str(s) for s in perturbation_function.free_symbols]

    configuration = PerturbationConfiguration(
        material=material,
        perturbation_function_holder=PerturbationFunctionHolder(
            function_str=perturbation_function,
            variables=variables,
        ),
        use_cartesian_coordinates=use_cartesian_coordinates,
    )
    if is_isometric:
        # builder = CellMatchingDistancePreservingSlabPerturbationBuilder()
        pass
    else:
        builder = PerturbationBuilder()
    return builder.get_material(configuration)
