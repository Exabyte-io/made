from typing import Any, List, Optional

import numpy as np
import sympy as sp
from mat3ra.made.material import Material
from mat3ra.made.periodic_table import get_atomic_mass_from_element
from pydantic import Field

from .elemental_function_holder import AtomicMassDependentFunctionHolder

BOLTZMANN_CONSTANT_EV_PER_K = 8.617333262145e-5


class MaxwellBoltzmannDisplacementHolder(AtomicMassDependentFunctionHolder):
    disorder_parameter: float = Field(exclude=True)
    kT: float = Field(exclude=True)
    random_state: Any = Field(default=None, exclude=True)
    is_mass_used: bool = Field(default=True, exclude=True)

    def __init__(
        self,
        atom_masses: List[float],
        disorder_parameter: float,
        random_seed: Optional[int] = None,
        is_mass_used: bool = True,
    ):
        if random_seed is not None:
            np.random.seed(random_seed)

        kT = BOLTZMANN_CONSTANT_EV_PER_K * disorder_parameter
        random_state = np.random.RandomState(random_seed) if random_seed is not None else np.random

        function_expr = sp.Symbol("f")
        super().__init__(
            function=function_expr,
            atom_masses=atom_masses,
            disorder_parameter=disorder_parameter,
            kT=kT,
            random_state=random_state,
            is_mass_used=is_mass_used,
        )

    def apply_function(self, coordinate: List[float], atom_index: int) -> List[float]:
        if atom_index < 0 or atom_index >= len(self.atom_masses):
            raise ValueError(f"Atom index {atom_index} out of range [0, {len(self.atom_masses)})")

        if self.is_mass_used:
            mass = self.atom_masses[atom_index]
            variance = self.kT / mass
        else:
            variance = self.kT

        std_dev = np.sqrt(variance)
        displacement = self.random_state.normal(0.0, std_dev, size=3)
        return displacement.tolist()


def create_maxwell_displacement_function(
    material: Material,
    disorder_parameter: float,
    random_seed: Optional[int] = None,
    is_mass_used: bool = True,
) -> MaxwellBoltzmannDisplacementHolder:
    """
    Create a Maxwell-Boltzmann displacement function for thermal perturbations.

    The function generates random 3D displacement vectors where each component
    follows a normal distribution with variance proportional to kT/m (if is_mass_used=True)
    or kT (if is_mass_used=False), where k is Boltzmann's constant, T is the disorder
    parameter, and m is atomic mass.

    Args:
        material: The material containing atoms to be perturbed.
        disorder_parameter: Disorder parameter (typically temperature in Kelvin for
                            Maxwell-Boltzmann distribution).
        random_seed: Optional random seed for deterministic behavior.
        is_mass_used: If True, displacement variance is kT/m (mass-dependent).
                     If False, displacement variance is kT (mass-independent).

    Returns:
        MaxwellBoltzmannDisplacementHolder that generates Maxwell-Boltzmann displacements.
    """
    atom_masses = []
    if is_mass_used:
        for element_value in material.basis.elements.values:
            atomic_mass = get_atomic_mass_from_element(element_value)
            atom_masses.append(atomic_mass)
    else:
        atom_masses = [1.0] * len(material.basis.elements.values)

    return MaxwellBoltzmannDisplacementHolder(
        atom_masses=atom_masses,
        disorder_parameter=disorder_parameter,
        random_seed=random_seed,
        is_mass_used=is_mass_used,
    )

