from typing import Union

import sympy as sp
from mat3ra.made.material import Material

from ..... import MaterialWithBuildMetadata
from .builders.base import PerturbationBuilder
from .builders.isometric import IsometricPerturbationBuilder
from .configuration import PerturbationConfiguration
from .functions import PerturbationFunctionHolder


def create_perturbation(
    material: Union[Material, MaterialWithBuildMetadata],
    perturbation_function: Union[sp.Expr, str],
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

    configuration = PerturbationConfiguration(
        material=material,
        perturbation_function_holder=PerturbationFunctionHolder(function=perturbation_function),
        use_cartesian_coordinates=use_cartesian_coordinates,
    )
    if is_isometric:
        builder = IsometricPerturbationBuilder()
    else:
        builder = PerturbationBuilder()
    return builder.get_material(configuration)
