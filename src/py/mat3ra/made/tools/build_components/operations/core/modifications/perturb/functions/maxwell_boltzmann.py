from typing import Any, Optional

import numpy as np
import sympy as sp
from mat3ra.made.periodic_table import get_atomic_mass_from_element
from pydantic import Field

from .elemental_function_holder import AtomicMassDependentFunctionHolder


class MaxwellBoltzmannDisplacementHolder(AtomicMassDependentFunctionHolder):
    disorder_parameter: float = Field(exclude=True)
    random_state: Any = Field(default=None, exclude=True)
    is_mass_used: bool = Field(default=True, exclude=True)

    def __init__(
        self,
        disorder_parameter: float,
        random_seed: Optional[int] = None,
        is_mass_used: bool = True,
    ):
        if random_seed is not None:
            np.random.seed(random_seed)

        random_state = np.random.RandomState(random_seed) if random_seed is not None else np.random

        function_expr = sp.Symbol("f")
        super().__init__(
            function=function_expr,
            disorder_parameter=disorder_parameter,
            random_state=random_state,
            is_mass_used=is_mass_used,
        )

    def apply_function(self, coordinate, material=None, atom_index: Optional[int] = None, **kwargs) -> list:
        if material is None or atom_index is None:
            raise ValueError("MaxwellBoltzmannDisplacementHolder requires 'material' and 'atom_index' kwargs")

        if self.is_mass_used:
            element = material.basis.elements.values[atom_index]
            mass = get_atomic_mass_from_element(element)
            variance = self.disorder_parameter / mass
        else:
            variance = self.disorder_parameter

        std_dev = np.sqrt(variance)
        displacement = self.random_state.normal(0.0, std_dev, size=3)
        return displacement.tolist()
