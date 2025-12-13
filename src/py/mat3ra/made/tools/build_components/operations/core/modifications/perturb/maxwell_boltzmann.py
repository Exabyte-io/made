from typing import Any, List, Optional

import numpy as np
import sympy as sp
from mat3ra.made.material import Material
from mat3ra.made.periodic_table import get_atomic_mass_from_element
from pydantic import Field

from .functions import ElementalFunctionHolder

BOLTZMANN_CONSTANT_EV_PER_K = 8.617333262145e-5


class MaxwellBoltzmannDisplacementHolder(ElementalFunctionHolder):
    temperature: float = Field(exclude=True)
    kT: float = Field(exclude=True)
    random_state: Any = Field(default=None, exclude=True)

    def __init__(
        self,
        atom_masses: List[float],
        temperature_in_kelvin: float,
        random_seed: Optional[int] = None,
    ):
        if random_seed is not None:
            np.random.seed(random_seed)

        temperature = temperature_in_kelvin
        kT = BOLTZMANN_CONSTANT_EV_PER_K * temperature_in_kelvin
        random_state = np.random.RandomState(random_seed) if random_seed is not None else np.random

        function_expr = sp.Symbol("f")
        super().__init__(
            function=function_expr,
            atom_masses=atom_masses,
            temperature=temperature,
            kT=kT,
            random_state=random_state,
        )

    def apply_function(self, coordinate: List[float], atom_index: int) -> List[float]:
        if atom_index < 0 or atom_index >= len(self.atom_masses):
            raise ValueError(f"Atom index {atom_index} out of range [0, {len(self.atom_masses)})")

        mass = self.atom_masses[atom_index]
        variance = self.kT / mass
        std_dev = np.sqrt(variance)
        displacement = self.random_state.normal(0.0, std_dev, size=3)
        return displacement.tolist()


def create_maxwell_displacement_function(
    material: Material, temperature_in_kelvin: float, random_seed: Optional[int] = None
) -> MaxwellBoltzmannDisplacementHolder:
    """
    Create a Maxwell-Boltzmann displacement function for thermal perturbations.

    The function generates random 3D displacement vectors where each component
    follows a normal distribution with variance proportional to kT/m, where
    k is Boltzmann's constant, T is temperature, and m is atomic mass.

    Args:
        material: The material containing atoms to be perturbed.
        temperature_in_kelvin: Temperature in Kelvin.
        random_seed: Optional random seed for deterministic behavior.

    Returns:
        MaxwellBoltzmannDisplacementHolder that generates Maxwell-Boltzmann displacements.
    """
    atom_masses = []
    for element_value in material.basis.elements.values:
        atomic_mass = get_atomic_mass_from_element(element_value)
        atom_masses.append(atomic_mass)

    return MaxwellBoltzmannDisplacementHolder(
        atom_masses=atom_masses,
        temperature_in_kelvin=temperature_in_kelvin,
        random_seed=random_seed,
    )

